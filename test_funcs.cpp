#ifdef _DEBUG
#include <iostream>
#include <vector>

#include "TimeCurve.h"
#include "MemoryManager.h"

#include "Mesh_1D2.h"
#include "Particle_1D_Mechanics.h"
#include "BoundaryCondition.h"
#include "OutputRequest.h"
#include "FileBuffer.h"

#include "TmpDataToHdf5.h"

void test_TimeCurve(void)
{
	size_t i;
	int div_num;
	double dt;
	double t, x, dx_dt;
	TimeCurveCubicSmooth_inline curve;
	curve.x_begin = 0.0;
	curve.x_end = 5.0;
	curve.t_len = 2.5;
	curve.init();

	div_num = 20;
	dt = curve.t_len / div_num;
	for (i = 0; i < (div_num + 1); i++)
	{
		t = i * dt;
		x = curve.x(t);
		dx_dt = curve.dx_dt(t);
		//std::cout << "t: " << t << "; x: " << x << "; dx_dt: " << dx_dt << std::endl;
	}
}

void test_FixedSizeMemoryManager(void)
{
	FixedSizeMemeory mem(sizeof(double));
	size_t i;
	double *tmp, *head, *iter;

	mem.reserve(sizeof(double)/2);
	for (i = 0; i < 20; i++)
	{
		tmp = (double *)mem.alloc();
		*tmp = (double)i;
	}
	//mem.clear();
	mem.compress();

	std::cout << "item num: " << mem.get_item_num() << std::endl;
	std::cout << "item size: " << mem.get_item_size() << std::endl;

	head = (double *)mem.get_mem();
	for (i = 0; i < mem.get_item_num(); i++)
		std::cout << *(double*)(mem.get_item_by_id(i)) << " ";
	//std::cout << head[i] << " ";
	std::cout << std::endl;

	for (iter = (double *)mem.get_first(); iter; iter = (double *)mem.get_next(iter))
		std::cout << *iter << " ";
	std::cout << std::endl;
}

void test_FlexSizeMemoryManager(void)
{
	FlexibleSizeMemory mem;
	size_t i;
	double *tmp, *iter;

	mem.reserve(1 * sizeof(double));
	for (i = 0; i < 20; i++)
	{
		tmp = (double *)mem.alloc(sizeof(double));
		*tmp = (double)i;
	}
	//mem.clear();
	mem.compress();

	for (iter = (double *)mem.get_first(); iter; iter = (double *)mem.get_next(iter))
	{
		std::cout << "address: " << iter << "; con: " << *iter
			<< "; size: " << mem.get_item_size(iter) << std::endl;
	}
}

void test_vaBC(void)
{
	size_t i;
	VelocityBCManager vbcs;
	AccelerationBCManager abcs;

	// velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 1;
	vbc.dof = DegreeOfFreedom::x;
	vbcs.add_bc(&vbc);
	VelocityBCParam_Constant vbc2;
	vbc2.index = 2;
	vbc2.dof = DegreeOfFreedom::x;
	vbc2.v = 1.0;
	vbcs.add_bc(&vbc2);
	VelocityBCParam_Linear vbc3;
	vbc3.index = 3;
	vbc3.dof = DegreeOfFreedom::x;
	vbc3.v_begin = 0.0;
	vbc3.v_end = 5.0;
	vbc3.t_length = 2.5;
	vbcs.add_bc(&vbc3);
	VelocityBCParam_CubicSmooth vbc4;
	vbc4.index = 4;
	vbc4.dof = DegreeOfFreedom::x;
	vbc4.v_begin = 0.0;
	vbc4.v_end = 5.0;
	vbc4.t_length = 2.5;
	vbcs.add_bc(&vbc4);
	
	// acceleration boundary conditions
	AccelerationBCParam_Fixed abc;
	abc.index = 1;
	abc.dof = DegreeOfFreedom::x;
	abcs.add_bc(&abc);
	AccelerationBCParam_Constant abc2;
	abc2.index = 2;
	abc2.dof = DegreeOfFreedom::x;
	abc2.a = 1.0;
	abcs.add_bc(&abc2);
	AccelerationBCParam_Linear abc3;
	abc3.index = 3;
	abc3.dof = DegreeOfFreedom::x;
	abc3.a_begin = 0.0;
	abc3.a_end = 5.0;
	abc3.t_length = 2.5;
	abcs.add_bc(&abc3);
	AccelerationBCParam_CubicSmooth abc4;
	abc4.index = 4;
	abc4.dof = DegreeOfFreedom::y;
	abc4.a_begin = 0.0;
	abc4.a_end = 5.0;
	abc4.t_length = 2.5;
	abcs.add_bc(&abc4);

	double tt[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };
	double tt_len = 6;
	double res_tmp;

	VelocityBC *vbc_iter;
	for (vbc_iter = vbcs.get_first(); vbc_iter; vbc_iter = vbcs.get_next(vbc_iter))
	{
		std::cout << "vbc n_id: " << vbc_iter->index
			<< "; dof: " << (unsigned int)vbc_iter->dof
			<< "; type: " << (unsigned int)vbc_iter->getType() << std::endl;
		for (i = 0; i < tt_len; i++)
		{
			res_tmp = vbc_iter->v(tt[i]);
			std::cout << res_tmp << " ";
		}
		std::cout << std::endl;
		for (i = 0; i < tt_len; i++)
		{
			res_tmp = vbc_iter->a(tt[i]);
			std::cout << res_tmp << " ";
		}
		std::cout << std::endl;
	}
	AccelerationBC *abc_iter;
	for (abc_iter = abcs.get_first(); abc_iter; abc_iter = abcs.get_next(abc_iter))
	{
		std::cout << "abc n_id: " << abc_iter->index
			<< "; doc: " << (unsigned int)abc_iter->dof
			<< "; type: " << (unsigned int)abc_iter->getType() << std::endl;
		for (i = 0; i < tt_len; i++)
			std::cout << abc_iter->a(tt[i]) << " ";
		std::cout << std::endl;
	}
}

void test_forceBC(void)
{
	size_t i;
	MassForceBCManager mfbcs;
	SurfaceForceBCManager sfbcs;

	// Mass force boundary conditions
	MassForceBCParam_Constant mfbc2;
	mfbc2.index = 2;
	mfbc2.dof = DegreeOfFreedom::x;
	mfbc2.massForce = 1.0;
	mfbcs.add_bc(&mfbc2);
	MassForceBCParam_Linear mfbc3;
	mfbc3.index = 3;
	mfbc3.dof = DegreeOfFreedom::x;
	mfbc3.massForce_begin = 0.0;
	mfbc3.massForce_end = 5.0;
	mfbc3.t_length = 2.5;
	mfbcs.add_bc(&mfbc3);
	MassForceBCParam_CubicSmooth mfbc4;
	mfbc4.index = 4;
	mfbc4.dof = DegreeOfFreedom::x;
	mfbc4.massForce_begin = 0.0;
	mfbc4.massForce_end = 5.0;
	mfbc4.t_length = 2.5;
	mfbcs.add_bc(&mfbc4);

	// Surface force boundary conditions
	SurfaceForceBCParam_Constant sfbc2;
	sfbc2.index = 2;
	sfbc2.dof = DegreeOfFreedom::x;
	sfbc2.surfaceForce = 1.0;
	sfbcs.add_bc(&sfbc2);
	SurfaceForceBCParam_Linear sfbc3;
	sfbc3.index = 3;
	sfbc3.dof = DegreeOfFreedom::x;
	sfbc3.surfaceForce_begin = 0.0;
	sfbc3.surfaceForce_end = 5.0;
	sfbc3.t_length = 2.5;
	sfbcs.add_bc(&sfbc3);
	SurfaceForceBCParam_CubicSmooth sfbc4;
	sfbc4.index = 4;
	sfbc4.dof = DegreeOfFreedom::x;
	sfbc4.surfaceForce_begin = 0.0;
	sfbc4.surfaceForce_end = 5.0;
	sfbc4.t_length = 2.5;
	sfbcs.add_bc(&sfbc4);

	double tt[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };
	double tt_len = 6;
	double res_tmp;

	MassForceBC *mfbc_iter;
	for (mfbc_iter = mfbcs.get_first(); mfbc_iter; mfbc_iter = mfbcs.get_next(mfbc_iter))
	{
		std::cout << "mfbc n_id: " << mfbc_iter->index
			<< "; dof: " << (unsigned int)mfbc_iter->dof
			<< "; type: " << (unsigned int)mfbc_iter->getType() << std::endl;
		for (i = 0; i < tt_len; i++)
		{
			res_tmp = mfbc_iter->massForce(tt[i]);
			std::cout << res_tmp << " ";
		}
		std::cout << std::endl;
	}
	SurfaceForceBC *sfbc_iter;
	for (sfbc_iter = sfbcs.get_first(); sfbc_iter; sfbc_iter = sfbcs.get_next(sfbc_iter))
	{
		std::cout << "sfbc n_id: " << sfbc_iter->index
			<< "; doc: " << (unsigned int)sfbc_iter->dof
			<< "; type: " << (unsigned int)sfbc_iter->getType() << std::endl;
		for (i = 0; i < tt_len; i++)
			std::cout << sfbc_iter->surfaceForce(tt[i]) << " ";
		std::cout << std::endl;
	}
}


void test_FileBuffer(void)
{
	size_t i;
	TmpDataFile f;
	FileBuffer fb;
	// Attribute buffer
	AttributeBuffer ab;
	// Data buffer
	DataBuffer<double> db;
	//HTmpData meta_data1;
	HTmpData meta_data2;
	HTmpData meta_data3;

	char data[] = "meta_datatag";
	char attr1[] = { 1, 1, 1, 1, 1 };
	char attr2[] = { 2, 2, 2, 2, 2 };
	char attr3[] = { 1, 2, 3, 4, 5, 6, 7 };
	char attr4[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	double *dat;

	f.initFile("ddd", 1);
	/*
	//meta_data1 = f.addMetaData("Test File Buffer", data, strlen(data));
	meta_data1 = f.addMetaData("Test File Buffer");
	fb.initBuffer(10, &f, meta_data1);
	fb.writeData(attr1, sizeof(attr1));
	fb.writeData(attr2, sizeof(attr2));
	fb.writeData(attr1, sizeof(attr1));
	fb.writeData(attr1, sizeof(attr1));
	fb.writeData(attr1, sizeof(attr1));
	fb.writeData(attr3, sizeof(attr3));
	//fb.addAttribute(attr2, sizeof(attr2));
	fb.flushBuffer();
	*/
	
	meta_data2 = f.addMetaData("Test Attribute Buffer");
	ab.initBuffer(15, &f, meta_data2);
	ab.addAttribute(attr1, sizeof(attr1));
	ab.addAttribute(attr4, sizeof(attr4));
	ab.addAttribute(attr1, sizeof(attr1));
	ab.addAttribute(attr2, sizeof(attr2));
	ab.addAttribute(attr1, sizeof(attr1));
	ab.addAttribute(attr3, sizeof(attr3));
	ab.flushBuffer();

	meta_data3 = f.addMetaData("Test Data Buffer");
	db.initBuffer(3, 2, &f, meta_data3);
	dat = db.getDataBuffer();
	for (i = 0; i < 3; i++)
		dat[i] = (double)i;
	dat = db.getDataBuffer();
	for (i = 0; i < 3; i++)
		dat[i] = (double)(i + 3);
	dat = db.getDataBuffer();
	for (i = 0; i < 3; i++)
		dat[i] = (double)(i + 6);
	db.flushBuffer();
}


void test_OutputRequest(void)
{
	double time_step = 1.0;

	Mesh_1D2 mesh;
	double coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
	mesh.initMesh(coords, sizeof(coords) / sizeof(coords[0]));
	// velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 1;
	vbc.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc);
	mesh.finish_init();

	std::vector<ObjectByParticle_1D_Mechanics> objects;
	ParticleParam_1D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// particle 1
	pcl_param.x = 0.5;
	pcl_param.mass = 1.0;
	pcl_param.density = 2.0;
	pcl_param.momentum1 = 0.0;
	pcl_param.stress11 = 0.0;
	pcl_param.strain11 = 0.0;
	pcl_param.estrain11 = 0.0;
	pcl_param.pstrain11 = 0.0;
	le_param.E = 100.0;
	le_param.nu = 0.0;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 2
	pcl_param.x = 1.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 3
	pcl_param.x = 2.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 4
	pcl_param.x = 3.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 5
	pcl_param.x = 4.5;
	objects[0].addParticle(&pcl_param, &le_param);
	/*
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 0.01;
	mfp_param.index = 1;
	objects[0].addMassForceBC(&mfp_param);
	// particle 2
	mfp_param.index = 2;
	objects[0].addMassForceBC(&mfp_param);
	// particle 3
	mfp_param.index = 3;
	objects[0].addMassForceBC(&mfp_param);
	// particle 4
	mfp_param.index = 4;
	objects[0].addMassForceBC(&mfp_param);
	// particle 5
	mfp_param.index = 5;
	objects[0].addMassForceBC(&mfp_param);
	// Surface forces
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x;
	sfp_param1.surfaceForce = 1.0;
	sfp_param1.index = 5;
	SurfaceForceBCParam_Linear sfp_param2;
	sfp_param2.dof = DegreeOfFreedom::x;
	sfp_param2.surfaceForce_begin = 0.0;
	sfp_param2.surfaceForce_end = 1.0;
	sfp_param2.t_length = time_step / 3.0 * 2.0;
	sfp_param2.index = 5;
	objects[0].addSurfaceForceBC(&sfp_param2);
	*/
	for (size_t i = 0; i < objects.size(); i++)
		objects[i].finish_init();

	OutputRequest output_request;
	Output *out1;
	Output *out2;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	output_request.setFileName("test_output");
	output_request.setStepTime(time_step);

	pcl_id.push_back(1);
	pcl_id.push_back(2);
	//pcl_id.push_back(3);
	//pcl_id.push_back(4);
	//pcl_id.push_back(5);
	//pcl_id.push_back(6);
	//pcl_id.push_back(7);
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::x));
	//field_id.push_back(2);
	//field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::isInMesh));
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::velocity1));
	//field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::stress11));
	//field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::strain11));
	
	out1 = output_request.addOutput("test_output1", 10);
	out1->addObject(&objects[0], &field_id, &pcl_id);
	//out1->addObject(&objects[0], &field_id);
	
	out2 = output_request.addOutput("test_output2", 15);
	out2->addObject(&objects[0], &field_id);

	output_request.finish_init();
	output_request.outputMeshGeometry(&mesh);

	output_request.output(1, 0.01);
	output_request.output(2, 0.05, true);
	output_request.output(3, 0.08);

	output_request.complete();

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


void test_TmpDataToHdf5(void)
{
	TmpDataToHdf5 trans;

	trans.init("test_result");
	trans.generate();

}

#endif // _DEBUG