#include <iostream>
#include <cassert>

#include "OutputRequest.h"
#include "Solver.h"

#include "Mesh.h"
#include "Mesh_1D2.h"
#include "Mesh_R2D4.h"

/* ---------------------------------- FileDataKey ------------------------------- */
FileDataKey::KeyType FileDataKey::OutputData = FileDataKey::KeyType("OutputData");
FileDataKey::KeyType FileDataKey::OutputAttribute = FileDataKey::KeyType("OutputAttribute");
FileDataKey::KeyType FileDataKey::OutputStepIndex = FileDataKey::KeyType("OutputStepIndex");
FileDataKey::KeyType FileDataKey::OutputStepName = FileDataKey::KeyType("OutputStepName");
FileDataKey::KeyType FileDataKey::OutputIndex = FileDataKey::KeyType("OutputIndex");
FileDataKey::KeyType FileDataKey::OutputName = FileDataKey::KeyType("OutputName");
FileDataKey::KeyType FileDataKey::ObjectNumber = FileDataKey::KeyType("ObjectNumber");
FileDataKey::KeyType FileDataKey::ObjectIndex = FileDataKey::KeyType("ObjectIndex");
FileDataKey::KeyType FileDataKey::ObjectName = FileDataKey::KeyType("ObjectName");
FileDataKey::KeyType FileDataKey::SimulationType = FileDataKey::KeyType("SimulationType");
FileDataKey::KeyType FileDataKey::SimulationTypeName = FileDataKey::KeyType("SimulationTypeName");
FileDataKey::KeyType FileDataKey::FieldNumber = FileDataKey::KeyType("FieldNumber");
FileDataKey::KeyType FileDataKey::FieldType = FileDataKey::KeyType("FieldType");
FileDataKey::KeyType FileDataKey::FieldName = FileDataKey::KeyType("FieldName");
FileDataKey::KeyType FileDataKey::ParticleNumber = FileDataKey::KeyType("ParticleNumber");
FileDataKey::KeyType FileDataKey::ParticleIndex = FileDataKey::KeyType("ParticleIndex");

FileDataKey::KeyType FileDataKey::BackgroundMeshAttribute = FileDataKey::KeyType("BackgroundMeshAttribute");
FileDataKey::KeyType FileDataKey::MeshType = FileDataKey::KeyType("MeshType");
FileDataKey::KeyType FileDataKey::MeshTypeName = FileDataKey::KeyType("MeshTypeName");
FileDataKey::KeyType FileDataKey::XCoordNum = FileDataKey::KeyType("XCoordNum");
FileDataKey::KeyType FileDataKey::XCoord = FileDataKey::KeyType("XCoord");
FileDataKey::KeyType FileDataKey::YCoordNum = FileDataKey::KeyType("YCoordNum");
FileDataKey::KeyType FileDataKey::YCoord = FileDataKey::KeyType("YCoord");

/* -------------------- ObjectOutputManager -------------------- */
// *** Need to modify this function when new simulation is added
int ObjectOutputManager::add_object_output(ObjectByParticle *obj,
	std::vector<unsigned long long> *flds, std::vector<size_t> *pcls)
{
	void *ptr;
	ObjectOutput_1D_Mechanics *out1;
	ObjectOutput_1D_Hydromechanics *out2;
	ObjectOutput_2D_Mechanics *out3;
	
	if (!obj || !flds) return -1;

	switch (obj->getSimulationType())
	{
	case SimulationType::Mechanics_1D:
		ptr = alloc(sizeof(ObjectOutput_1D_Mechanics));
		if (!ptr) return -3;
		out1 = new (ptr) ObjectOutput_1D_Mechanics(
			static_cast<ObjectByParticle_1D_Mechanics*>(obj));
		out1->addFields(flds);
		out1->addParticles(pcls);
		break;
	case SimulationType::Hydromechanics_1D:
		ptr = alloc(sizeof(ObjectOutput_1D_Hydromechanics));
		if (!ptr) return -3;
		out2 = new (ptr) ObjectOutput_1D_Hydromechanics(
			static_cast<ObjectByParticle_1D_Hydromechanics*>(obj));
		out2->addFields(flds);
		out2->addParticles(pcls);
		break;
	case SimulationType::Mechanics_2D:
		ptr = alloc(sizeof(ObjectOutput_2D_Mechanics));
		if (!ptr) return -3;
		out3 = new (ptr) ObjectOutput_2D_Mechanics(
			static_cast<ObjectByParticle_2D_Mechanics*>(obj));
		out3->addFields(flds);
		out3->addParticles(pcls);
		break;
	default:
		return -2;
	}

	return 0;
}

ObjectOutputManager::~ObjectOutputManager()
{
	ObjectOutput *out_iter;
	for (out_iter = static_cast<ObjectOutput *>(get_first()); out_iter;
		 out_iter = static_cast<ObjectOutput *>(get_next(out_iter)))
		out_iter->~ObjectOutput();
}

/* --------------------------- Output --------------------------- */
size_t Output::curIndex = 0;

Output::Output(const char *out_name, size_t num, Output *nt, size_t bsrv) :
	index(++curIndex), output_num(num), next(nt),
	t_time(0.0), t_inc(0.0), t_tol(0.0),
	object_num(0),
	buffer_size_recomended_value(bsrv)
{
	if (out_name) name = out_name;
}

int Output::init(TmpDataFile &outFile, Solver &solver)
{
	size_t i;
	HTmpData attr_handle, data_handle;
	ObjectOutput *obj_iter;
	size_t row_len, row_num;
	char meta_name_tmp[64];
	const char *name_tmp;
	unsigned long long num_tmp;

	// Initialize time increment
	t_inc = solver.t_step / output_num;
	t_tol = t_inc * 0.01;
	t_time = t_inc - t_tol + solver.t_start;

	// Output attribute
	memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
	sprintf(meta_name_tmp, "%s%llu%s%llu", FileDataKey::OutputAttribute.key,
		(unsigned long long)index, FileDataKey::OutputStepIndex.key, solver.index);
	attr_handle = outFile.addMetaData(meta_name_tmp);
	attr_buffer.initBuffer(256, &outFile, attr_handle);
	// output index
	attr_buffer.addAttribute(FileDataKey::OutputIndex.key, FileDataKey::OutputIndex.len);
	num_tmp = (unsigned long long)index;
	attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
	// output name
	attr_buffer.addAttribute(FileDataKey::OutputName.key, FileDataKey::OutputName.len);
	if (name.size()) name_tmp = name.c_str();
	else
	{
		memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
		sprintf(meta_name_tmp, "Output%llu", num_tmp);
		name_tmp = meta_name_tmp;
	}
	attr_buffer.addAttribute(name_tmp, strlen(name_tmp));
	// step index
	attr_buffer.addAttribute(FileDataKey::OutputStepIndex.key, FileDataKey::OutputStepIndex.len);
	num_tmp = (size_t)solver.index;
	attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
	// step name
	attr_buffer.addAttribute(FileDataKey::OutputStepName.key, FileDataKey::OutputStepName.len);
	if (solver.name.size()) name_tmp = solver.name.c_str();
	else
	{
		memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
		sprintf(meta_name_tmp, "Step%llu", num_tmp);
		name_tmp = meta_name_tmp;
	}
	attr_buffer.addAttribute(name_tmp, strlen(name_tmp));

	// object num
	attr_buffer.addAttribute(FileDataKey::ObjectNumber.key, FileDataKey::ObjectNumber.len);
	num_tmp = (unsigned long long)object_num;
	attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
	// for each object
	for (obj_iter = output_mem.get_first(); obj_iter;
		 obj_iter = output_mem.get_next(obj_iter))
	{
		// object index
		attr_buffer.addAttribute(FileDataKey::ObjectIndex.key, FileDataKey::ObjectIndex.len);
		num_tmp = (unsigned long long)obj_iter->getObjectIndex();
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// object name
		attr_buffer.addAttribute(FileDataKey::ObjectName.key, FileDataKey::ObjectName.len);
		if (obj_iter->getObjectName().size())
			name_tmp = obj_iter->getObjectName().c_str();
		else
		{
			memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
			sprintf(meta_name_tmp, "Object%llu", num_tmp);
			name_tmp = meta_name_tmp;
		}
		attr_buffer.addAttribute(name_tmp, strlen(name_tmp));
		// simulation type
		attr_buffer.addAttribute(FileDataKey::SimulationType.key, FileDataKey::SimulationType.len);
		num_tmp = (unsigned long long)obj_iter->getSimulationType();
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// simulation type name
		attr_buffer.addAttribute(FileDataKey::SimulationTypeName.key, FileDataKey::SimulationTypeName.len);
		name_tmp = SimulationTypeName::getName((size_t)num_tmp);
		attr_buffer.addAttribute(name_tmp, strlen(name_tmp));
		// field num
		attr_buffer.addAttribute(FileDataKey::FieldNumber.key, FileDataKey::FieldNumber.len);
		num_tmp = (unsigned long long)(obj_iter->field_type_num);
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// output fields
		// field type
		attr_buffer.addAttribute(FileDataKey::FieldType.key, FileDataKey::FieldType.len);
		for (i = 0; i < obj_iter->field_type_num; i++)
		{
			num_tmp = (unsigned long long)(obj_iter->field_type_array[i]);
			attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		}
		// field name
		attr_buffer.addAttribute(FileDataKey::FieldName.key, FileDataKey::FieldName.len);
		for (i = 0; i < obj_iter->field_type_num; i++)
		{
			num_tmp = (unsigned long long)(obj_iter->field_type_array[i]);
			name_tmp = obj_iter->getFieldName((size_t)num_tmp);
			if (!name_tmp) // if no default field name
			{
				memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
				sprintf(meta_name_tmp, "Field%llu", num_tmp);
				name_tmp = meta_name_tmp;
			}
			attr_buffer.addAttribute(name_tmp, strlen(name_tmp));
		}
		// particles num
		attr_buffer.addAttribute(FileDataKey::ParticleNumber.key, FileDataKey::ParticleNumber.len);
		num_tmp = (unsigned long long)(obj_iter->particle_index_num);
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// output particles
		// particle index
		attr_buffer.addAttribute(FileDataKey::ParticleIndex.key, FileDataKey::ParticleIndex.len);
		for (i = 0; i < obj_iter->particle_index_num; i++)
		{
			num_tmp = (unsigned long long)(obj_iter->particle_index_array[i]);
			attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		}
	}
	// Finish attribute output
	attr_buffer.flushBuffer();

	// Initialize buffer for data output
	// calculate the total number of output and initialize each ObjectOutput class
	row_len = 2; // for interation index and time
	for (obj_iter = output_mem.get_first(); obj_iter;
		obj_iter = output_mem.get_next(obj_iter))
	{
		obj_iter->init();
		row_len += obj_iter->getOutputNum();
	}
	row_num = buffer_size_recomended_value / (row_len * sizeof(double));
	if (!row_num) row_num = 1;
	memset(meta_name_tmp, 0, sizeof(meta_name_tmp));
	sprintf(meta_name_tmp, "%s%llu%s%llu", FileDataKey::OutputData.key,
		(unsigned long long)index, FileDataKey::OutputStepIndex.key, solver.index);
	data_handle = outFile.addMetaData(meta_name_tmp);
	data_buffer.initBuffer(row_len, row_num, &outFile, data_handle);

	return 0;
}

#define OUTPUT_DATA_CODE(iter_index, t_cur)            \
	do {                                               \
		double *data_row_buffer;                       \
		ObjectOutput *iter;                            \
		union                                          \
		{                                              \
			unsigned long long ull;                    \
			double d;                                  \
		} ulltod_tmp;                                  \
		data_row_buffer = data_buffer.getDataBuffer(); \
		/* output iteraction index */                  \
		ulltod_tmp.ull = iter_index;                   \
		data_row_buffer[0] = ulltod_tmp.d;             \
		/* output time */                              \
		data_row_buffer[1] = t_cur;                    \
		data_row_buffer += 2;                          \
		/* output data of each objects */              \
		for (iter = output_mem.get_first(); iter; iter = output_mem.get_next(iter)) \
			data_row_buffer = iter->output(data_row_buffer); \
	} while(0)

void Output::output_f(unsigned long long iter_index, double t_cur)
{
	OUTPUT_DATA_CODE(iter_index, t_cur);
}

double Output::output(unsigned long long iter_index, double t_cur)
{
	if (t_cur < t_time) return t_time;

	// update time for next output
	t_time += t_inc;
	if (t_time < t_cur) t_time = t_cur + t_inc - t_tol;

	OUTPUT_DATA_CODE(iter_index, t_cur);

	return t_time;
}

void Output::complete() { data_buffer.flushBuffer(); }


/* ----------------------- OutputRequest ----------------------- */
OutputRequest::OutputRequest(const char *fname) :
	file_name(fname),
	next_output_time(0.0),
	top(nullptr), tail(nullptr), output_num(0)
{
	if (file_name.length())
	{
		if (file.initFile(file_name.c_str(), 1))
			throw std::exception("Fail to init file of OutputRequest.");
		file_name = file.getFileName();
	}
}

Output *OutputRequest::addOutput(const char *out_name,
	size_t output_num, size_t buf_size)
{
	Output *tmp;

	tmp = new Output(out_name, output_num, nullptr, buf_size);
	if (!tmp) return nullptr;
	if (tail) tail->next = tmp;
	else top = tmp, tail = tmp;
	++output_num;

	return tmp;
}


int OutputRequest::deleteOutput(const char *out_name)
{
	Output *tmp, *prev_tmp;

	tmp = top;
	prev_tmp = nullptr;
	while (tmp)
	{
		if (!strcmp(tmp->name.c_str(), out_name))
		{
			if (prev_tmp) // tmp is not the first node
			{
				prev_tmp->next = tmp->next;
				delete tmp;
			}
			else
			{
				top->next = tmp->next;
				delete tmp;
			}
			--output_num;
			return 0;
		}
		prev_tmp = tmp;
		tmp = top->next;
	}
	
	return -1;
}


void OutputRequest::output_data_header(void)
{
	Output *tmp;
	size_t step_index = solver->index;
	double step_time = solver->t_step;
	double step_start_time = solver->t_start;

	tmp = top;
	while (tmp)
	{
		tmp->init(file, *solver);
		tmp = tmp->next;
	}
}


// Force to output regardless of the time
void OutputRequest::output_f(void)
{
	Output *tmp;
	unsigned long long iter_index = solver->iteration_index;
	double cur_time = solver->t_cur;

	tmp = top;
	while (tmp)
	{
		tmp->output_f(iter_index, cur_time);
		tmp = tmp->next;
	}

	return;
}

void OutputRequest::output(void)
{
	Output *tmp;
	double time_tmp;
	unsigned long long iter_index;
	double cur_time = solver->t_cur;

	if (cur_time > next_output_time)
	{
		if (!top) return;
		iter_index = solver->iteration_index;
		next_output_time = top->output(iter_index, cur_time);
		tmp = top->next;
		while (tmp)
		{
			time_tmp = tmp->output(iter_index, cur_time);
			next_output_time = time_tmp < next_output_time ? time_tmp : next_output_time;
			tmp = tmp->next;
		}
	}
}

void OutputRequest::complete()
{
	Output *tmp;
	tmp = top;
	while (tmp)
	{
		tmp->complete();
		tmp = tmp->next;
	}
	file.flushFile();
}

OutputRequest::~OutputRequest() { reset(); }

void OutputRequest::reset()
{
	Output *tmp;
	while (top)
	{
		tmp = top;
		top = top->next;
		delete tmp;
	}
}

void OutputRequest::outputMeshGeometry(Mesh *mh)
{
	size_t i;
	double *nodeCoords;
	Mesh_1D2 * mh1;
	Mesh_R2D4 *mh2;
	AttributeBuffer attr_buffer;
	HTmpData attr_handle;
	unsigned long long num_tmp;
	const char *name_tmp;

	attr_handle = file .addMetaData(FileDataKey::BackgroundMeshAttribute.key);
	attr_buffer.initBuffer(256, &file, attr_handle);
	// mesh type index
	attr_buffer.addAttribute(FileDataKey::MeshType.key, FileDataKey::MeshType.len);
	num_tmp = (unsigned long long)mh->getMeshType();
	attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
	// mesh type name
	attr_buffer.addAttribute(FileDataKey::MeshTypeName.key, FileDataKey::MeshTypeName.len);
	name_tmp = MeshTypeName::getName((size_t)num_tmp);
	attr_buffer.addAttribute(name_tmp, strlen(name_tmp));

	// Output coordinates
	switch (mh->getMeshType())
	{
	case MeshType::Mesh_1D2:
		mh1 = static_cast<Mesh_1D2 *>(mh);
		// number of coordinates
		attr_buffer.addAttribute(FileDataKey::XCoordNum.key, FileDataKey::XCoordNum.len);
		num_tmp = (unsigned long long)mh1->getNodeNum();
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// coorindates
		attr_buffer.addAttribute(FileDataKey::XCoord.key, FileDataKey::XCoord.len);
		nodeCoords = (double *)mh1->getNodeCoords();
		for (i = 0; i < num_tmp; i++)
			attr_buffer.addAttribute(nodeCoords + i, sizeof(double));
		break;
	case MeshType::Mesh_R2D4:
		mh2 = static_cast<Mesh_R2D4 *>(mh);
		// number of x coordinates
		attr_buffer.addAttribute(FileDataKey::XCoordNum.key, FileDataKey::XCoordNum.len);
		num_tmp = (unsigned long long)mh2->getXCoordNum();
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// x coordinates
		attr_buffer.addAttribute(FileDataKey::XCoord.key, FileDataKey::XCoord.len);
		nodeCoords = (double *)mh2->getXNodeCoords();
		for (i = 0; i < num_tmp; i++)
			attr_buffer.addAttribute(nodeCoords + i, sizeof(double));
		// number of y coordinates
		attr_buffer.addAttribute(FileDataKey::YCoordNum.key, FileDataKey::YCoordNum.len);
		num_tmp = (unsigned long long)mh2->getYCoordNum();
		attr_buffer.addAttribute(&num_tmp, sizeof(unsigned long long));
		// y coordinates
		attr_buffer.addAttribute(FileDataKey::YCoord.key, FileDataKey::YCoord.len);
		nodeCoords = (double *)mh2->getYNodeCoords();
		for (i = 0; i < num_tmp; i++)
			attr_buffer.addAttribute(nodeCoords + i, sizeof(double));
		break;
	default:
		break;
	}

	// finish output
	attr_buffer.flushBuffer();
}

