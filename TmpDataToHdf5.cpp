#include <iostream>

#include "OutputRequest.h"
#include "Mesh.h"
#include "TmpDataToHdf5.h"

TmpDataToHdf5::TmpDataToHdf5() : hdf5_file_id(-1) {}

TmpDataToHdf5::~TmpDataToHdf5() { complete(); }

void TmpDataToHdf5::complete(void)
{
	// close TmpData file
	file.closeFile();

	/*// close all output group
	std::map<unsigned long long, hid_t>::iterator out_map_iter = output_id_map.begin();
	while (out_map_iter != output_id_map.end())
	{
		H5Gclose(out_map_iter->second);
		++out_map_iter;
	}
	output_id_map.clear();*/

	// close hdf5 file
	if (hdf5_file_id >= 0)
	{
		H5Fclose(hdf5_file_id);
		hdf5_file_id = -1;
	}
}

std::string formHdf5Name(std::string &file_name)
{
	char *name_tmp;
	size_t name_len_tmp;
	std::string res;

	name_len_tmp = file_name.size();
	name_tmp = new char[name_len_tmp + 1];
	if (!name_tmp) return std::string();
	strcpy(name_tmp, file_name.c_str());
	strcpy(name_tmp + name_len_tmp - 8, ".h5");
	res = name_tmp;
	delete name_tmp;
	return res;
}

int TmpDataToHdf5::init(const char *file_name)
{
	int res;

	// Open TmpData file
	res = file.initFile(file_name, 0);
	if (res) return res;
	// TmpData file name
	tmp_data_file_name = file.getFileName();
	// hdf5 file name
	hdf5_file_name = formHdf5Name(tmp_data_file_name);
	//std::cout << hdf5_file_name << std::endl;
	// Create hdf5 file
	hdf5_file_id = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (hdf5_file_id < 0) return -1;

	return 0;
}

int TmpDataToHdf5::init(OutputRequest &out_req)
{
	int res;
	
	out_req.file.closeFile();
	// TmpData file name
	tmp_data_file_name = out_req.file.getFileName();
	// Open TmpData file
	res = file.initFile(tmp_data_file_name.c_str(), 0);
	if (res) return res;
	
	// hdf5 file name
	hdf5_file_name = formHdf5Name(tmp_data_file_name);
	//std::cout << hdf5_file_name << std::endl;
	// Create hdf5 file
	hdf5_file_id = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (hdf5_file_id < 0) return -1;
	
	return 0;
}

int TmpDataToHdf5::generate(void)
{

	//hid_t out_grp_id; // group id for output
	//hid_t obj_grp_id; // group id for object
	//hid_t field_id;   // dataset id for field
	hid_t datatype_id;
	//hid_t mem_dataspace_id, file_dataspace_id;
	//hid_t property_id;
	//hid_t attribute_id;
	//hid_t attribute_dataspace_id;
	//size_t attribute_len_tmp; // used to init attributes

	HTmpData meta_handle;
	const char *data_type;


	// Initialize loop
	// set type of the field data
	datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
	H5Tset_order(datatype_id, H5T_ORDER_LE);
	while (true) // for each field (metadata)
	{
		// read each meta data of TmpData File
		meta_handle = file.getNextMetaData();
		if (!meta_handle) break;

		/*
			Note that currently the output data should only appear after
		its corresponding attribute.
		*/
		data_type = file.getDataType(meta_handle);
		//std::cout << data_type << std::endl;
		if (!strncmp(data_type, FileDataKey::OutputAttribute.key, FileDataKey::OutputAttribute.len))
		{
			// Output attribute
			parse_output_attribute(meta_handle);
		}
		else if (!strncmp(data_type, FileDataKey::OutputData.key, FileDataKey::OutputData.len))
		{
			// Output data
			parse_output_data(meta_handle);
		}
		else if (!strncmp(data_type, FileDataKey::BackgroundMeshAttribute.key, FileDataKey::BackgroundMeshAttribute.len))
		{
			// Background mesh attribute
			parse_background_mesh_attribute(meta_handle);
		}
	}


	return 0;
}


class Token
{
protected:
	char *str;
	size_t len;
	size_t capacity;
public:
	Token() : str(nullptr), len(0), capacity(0) {}
	Token(const void *ex_str, size_t ex_str_len) { set_buf(ex_str, ex_str_len); }
	Token(const char *ex_str) { set_buf(ex_str, strlen(ex_str)); }
	~Token() { if (str) free(str); }
	void *resize(size_t new_size)
	{
		char *tmp;
		if (new_size > capacity)
		{
			tmp = (char *)realloc(str, new_size);
			if (tmp) // alloc succeed
			{
				capacity = new_size;
				str = tmp;
			}
			else // realloc fail
			{
				capacity = 0;
				if (str) free(str);
				str = nullptr;
			}
		}
		return str;
	}
	inline size_t get_len(void) noexcept { return len; }
	inline void *get_buf(void) noexcept { return str; }
	void *set_buf(const void* ex_str, size_t ex_str_len)
	{
		if (ex_str_len + 1 > capacity)
			if (!resize(ex_str_len + 1)) return nullptr;
		memcpy(str, ex_str, ex_str_len);
		str[ex_str_len] = '\0';
		len = ex_str_len;
		return str;
	}
	void clear()
	{
		if (str) free(str);
		str = nullptr;
		len = 0;
		capacity = 0;
	}
};

const void *getToken(const void *data, Token &tk)
{
	unsigned long long data_len;
	data_len = *(unsigned long long*)data;
	if (tk.set_buf((char *)data + sizeof(unsigned long long), (size_t)data_len))
		return (char *)data + sizeof(unsigned long long) + (size_t)data_len;
	return nullptr;
}

void add_attr(hid_t target, const char *attr_name, unsigned long long attr)
{
	hid_t attribute_dataspace_id, attribute_id;
	attribute_dataspace_id = H5Screate(H5S_SCALAR);
	attribute_id = H5Acreate(target, attr_name, H5T_NATIVE_ULLONG,
		attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &attr);
	H5Aclose(attribute_id);
	H5Sclose(attribute_dataspace_id);
}

void add_attr(hid_t target, const char *attr_name, const char *attr)
{
	hid_t attribute_dataspace_id, attribute_id;
	hsize_t attribute_len_tmp = strlen(attr);
	if (attribute_len_tmp)
	{
		attribute_dataspace_id = H5Screate_simple(1, &attribute_len_tmp, NULL);
		attribute_id = H5Acreate(target, attr_name, H5T_NATIVE_CHAR,
			attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(attribute_id, H5T_NATIVE_CHAR, attr);
		H5Aclose(attribute_id);
		H5Sclose(attribute_dataspace_id);
	}
}

int TmpDataToHdf5::parse_output_attribute(HTmpData meta_handle)
{
	size_t i, j;
	Token data_buffer_tmp;
	size_t data_len;
	char *data_buffer;
	Token tk_tmp;

	ObjectInfo *cur_obj_info;

	hid_t out_grp_id;
	hid_t obj_grp_id;
	hid_t field_dataset_id;
	hid_t field_datatype_id;
	hid_t field_dataspace_id;
	hsize_t dims[2], max_dims[2];
	hid_t property_id;
	hsize_t chunk_size;
	hsize_t chunk_dims[2];

	unsigned long long out_index, step_index, obj_num;
	unsigned long long obj_index, sim_type, fld_num, pcl_num;
	unsigned long long fld_type, pcl_index;

	hid_t pcl_datatype_id;
	hid_t pcl_dataspace_id;
	hid_t pcl_dataset_id;
	std::vector<unsigned long long> pcl_index_array;
	hsize_t pcl_num_hsize;

	std::string key_name, key_name2, key_name3, key_name4;
	std::string output_name, step_name;
	std::string output_grp_name;
	std::string key_fld_type, key_fld_name;

	// Assume all attribute have been written in only one chunk
	// get data
	data_len = file.getNextDataChunk(meta_handle);
	if (!data_len) return -1;
	data_buffer_tmp.resize(data_len);
	//std::cout << "data_buffer len: " << data_len << std::endl;
	data_buffer = (char *)data_buffer_tmp.get_buf();
	file.getData(meta_handle, data_buffer);
	
	tk_tmp.resize(128);

	// output index
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	out_index = *(unsigned long long *)(tk_tmp.get_buf());
	//std::cout << out_index << std::endl;

	// output name
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name2 = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	output_name = (char *)tk_tmp.get_buf();
	//std::cout << (char *)out_name.get_buf() << std::endl;

	// step index
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name3 = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	step_index = *(unsigned long long *)(tk_tmp.get_buf());
	//std::cout << step_index << std::endl;

	// step name
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name4 = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	step_name = (char *)tk_tmp.get_buf();
	//std::cout << step_name << std::endl;

	// create group for this output
	output_grp_name = output_name + std::string("-") + step_name;
	out_grp_id = H5Gcreate(hdf5_file_id, output_grp_name.c_str(),
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// write output index attribute
	add_attr(out_grp_id, key_name.c_str(), out_index);
	// write output name attribute
	add_attr(out_grp_id, key_name2.c_str(), output_name.c_str());
	// write step index attribute
	add_attr(out_grp_id, key_name3.c_str(), step_index);
	// write step name attribute
	add_attr(out_grp_id, key_name4.c_str(), step_name.c_str());

	cur_output_id = out_index;
	cur_output_name = output_grp_name.c_str();

	// object number
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	obj_num = *(unsigned long long *)tk_tmp.get_buf();
	//std::cout << obj_num << std::endl;
	add_attr(out_grp_id, key_name.c_str(), obj_num);

	cur_output_info.clear();
	// for each objects
	for (i = 0; i < obj_num; i++)
	{
		// object index
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		obj_index = *(unsigned long long *)tk_tmp.get_buf();
		//std::cout << "obj index: " << obj_index << std::endl;
		// object name
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name2 = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << "obj name: " << (char *)tk_tmp.get_buf() << std::endl;
		// Create group
		obj_grp_id = H5Gcreate(out_grp_id, (char *)tk_tmp.get_buf(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		// write object index attribute
		add_attr(obj_grp_id, key_name.c_str(), obj_index);
		// write object name attribute
		add_attr(obj_grp_id, key_name2.c_str(), (char *)tk_tmp.get_buf());

		cur_output_info.push_back(ObjectInfo());
		cur_obj_info = &cur_output_info.back();
		cur_obj_info->obj_index = obj_index;
		cur_obj_info->obj_name = (char *)tk_tmp.get_buf();

		// simulation type
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		sim_type = *(unsigned long long *)tk_tmp.get_buf();
		//std::cout << "sim type: " << sim_type << std::endl;
		add_attr(obj_grp_id, key_name.c_str(), sim_type);
		
		// simulation type name
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << "sim name: " << (char *)tk_tmp.get_buf() << std::endl;
		add_attr(obj_grp_id, key_name.c_str(), (char *)tk_tmp.get_buf());

		// field number
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		fld_num = *(unsigned long long *)tk_tmp.get_buf();
		//std::cout << "fld num: " << fld_num << std::endl;
		add_attr(obj_grp_id, key_name.c_str(), fld_num);

		cur_obj_info->fld_num = fld_num;

		// Read fields info
		// field type
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_fld_type = (char *)tk_tmp.get_buf();
		for (j = 0; j < fld_num; j++)
		{
			data_buffer = (char *)getToken(data_buffer, tk_tmp);
			fld_type = *(unsigned long long *)tk_tmp.get_buf();
			//std::cout << "fld type: " << fld_type << std::endl;

			cur_obj_info->fld_type.push_back(fld_type);
		}
		// field name
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_fld_name = (char *)tk_tmp.get_buf();
		for (j = 0; j < fld_num; j++)
		{
			data_buffer = (char *)getToken(data_buffer, tk_tmp);
			//std::cout << "fld name: " << (char *)tk_tmp.get_buf() << std::endl;
			
			cur_obj_info->fld_name.push_back(std::string((char *)tk_tmp.get_buf()));
		}

		// particle number
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		pcl_num = *(unsigned long long *)tk_tmp.get_buf();
		//std::cout << "pcl num: " << pcl_num << std::endl;
		add_attr(obj_grp_id, key_name.c_str(), pcl_num);

		cur_obj_info->pcl_num = pcl_num;

		// particle index
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		for (j = 0; j < pcl_num; j++)
		{
			data_buffer = (char *)getToken(data_buffer, tk_tmp);
			pcl_index = *(unsigned long long *)tk_tmp.get_buf();
			//std::cout << pcl_index << " ";
			pcl_index_array.push_back(pcl_index);
		}
		//std::cout << std::endl;
		// Create dataset for particle index
		if (pcl_num)
		{
			// 1. Create data type (double)
			pcl_datatype_id = H5Tcopy(H5T_NATIVE_ULLONG);
			H5Tset_order(pcl_datatype_id, H5T_ORDER_LE);
			// 2. Create dataspace
			pcl_num_hsize = (hsize_t)pcl_num;
			pcl_dataspace_id = H5Screate_simple(1, &pcl_num_hsize, NULL);
			// 3. Create dataset and write data
			pcl_dataset_id = H5Dcreate(obj_grp_id, key_name.c_str(),
				pcl_datatype_id, pcl_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Dwrite(pcl_dataset_id, pcl_datatype_id, pcl_dataspace_id, H5S_ALL, H5P_DEFAULT, &pcl_index_array[0]);
			pcl_index_array.clear();
			H5Dclose(pcl_dataset_id);
			// 4. Clear
			H5Tclose(pcl_datatype_id);
			H5Sclose(pcl_dataspace_id);
		}

		// Add dataset to obj_grp
		// 1. Create data type (double)
		field_datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
		H5Tset_order(field_datatype_id, H5T_ORDER_LE);
		// 2. Create dataspace
		dims[0] = 0; // no data initially
		dims[1] = pcl_num;
		max_dims[0] = H5S_UNLIMITED;
		max_dims[1] = pcl_num;
		field_dataspace_id = H5Screate_simple(2, dims, max_dims);
		// 3. Set dataset proporty
		property_id = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_layout(property_id, H5D_CHUNKED);
		// Assume each chunk is around 50kB
		chunk_size = 50 * 1024 / sizeof(double) / pcl_num;
		chunk_size = chunk_size ? chunk_size : 1; // chunk size > 1
		chunk_dims[0] = chunk_size;
		chunk_dims[1] = dims[1];
		H5Pset_chunk(property_id, 2, chunk_dims);
		// 4. Create dataset and add attributes
		for (j = 0; j < fld_num; j++)
		{
			// Create dataset
			field_dataset_id = H5Dcreate(obj_grp_id, cur_obj_info->fld_name[j].c_str(), 
				field_datatype_id, field_dataspace_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
			
			// Add field attributes
			// field type
			add_attr(field_dataset_id, key_fld_type.c_str(), cur_obj_info->fld_type[j]);
			// field type name
			add_attr(field_dataset_id, key_fld_name.c_str(), cur_obj_info->fld_name[j].c_str());
			
			H5Dclose(field_dataset_id);
		}
		// 5. Clear
		H5Tclose(field_datatype_id);
		H5Sclose(field_dataspace_id);
		H5Pclose(property_id);

		H5Gclose(obj_grp_id);
	}

	H5Gclose(out_grp_id);

	return 0;
}


int TmpDataToHdf5::parse_output_data(HTmpData meta_handle)
{
	size_t i, j;
	Token data_buffer_tmp;
	size_t data_len;
	double *data_buffer; // used to read data chunk
	ObjectInfo *cur_obj_info;

	hid_t out_grp_id;
	hid_t obj_grp_id;
	hid_t datatype_id;
	hid_t fld_dataset_id;
	hid_t mem_dataspace_id;
	hid_t file_dataspace_id;
	hsize_t dims[2], start[2], count[2];
	hid_t property_id;

	hid_t time_dataset_id;
	hid_t time_dataspace_id;
	hsize_t time_num, max_time_num, time_chunk_size;
	hsize_t time_start, time_count;

	hid_t iter_index_dataset_id;
	hid_t ull_datatype_id;
	hid_t iter_index_dataspace_id;
	hsize_t iter_index_num, max_iter_index_num, iter_index_chunk_size;
	hsize_t iter_index_start, iter_index_count;

	unsigned long long obj_num, fld_num, pcl_num;

	size_t total_column_num;
	size_t total_row_num, chunk_row_num;

	Token tmp_data_buf;
	double *tmp_data;
	size_t data_index_i, data_index_j; // used for iteration
	size_t cur_column;

	datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
	H5Tset_order(datatype_id, H5T_ORDER_LE);

	//std::cout << cur_output_name.c_str() << std::endl;
	out_grp_id = H5Gopen(hdf5_file_id, cur_output_name.c_str(), H5P_DEFAULT);
	
	// Create dataset for iteration index
	// Create datatype for unsigned long long
	ull_datatype_id = H5Tcopy(H5T_NATIVE_ULLONG);
	H5Tset_order(ull_datatype_id, H5T_ORDER_LE);
	// Create dataspace
	iter_index_num = 0;
	max_iter_index_num = H5S_UNLIMITED;
	iter_index_dataspace_id = H5Screate_simple(1, &iter_index_num, &max_iter_index_num);
	// Set dataset proporty
	property_id = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(property_id, H5D_CHUNKED);
	// Assume each chunk is around 0.5K
	iter_index_chunk_size = 512 / sizeof(unsigned long long);
	H5Pset_chunk(property_id, 1, &iter_index_chunk_size);
	// Create dataset
	iter_index_dataset_id = H5Dcreate(out_grp_id, "IterationIndex", 
		ull_datatype_id, iter_index_dataspace_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
	// Clean
	H5Sclose(iter_index_dataspace_id);
	H5Pclose(property_id);

	// Create dataset for time
	// Create dataspace
	time_num = 0;
	max_time_num = H5S_UNLIMITED;
	time_dataspace_id = H5Screate_simple(1, &time_num, &max_time_num);
	// Set dataset proporty
	property_id = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(property_id, H5D_CHUNKED);
	// Assume each chunk is around 0.5K
	time_chunk_size = 512 / sizeof(double);
	H5Pset_chunk(property_id, 1, &time_chunk_size);
	// Create dataset
	time_dataset_id = H5Dcreate(out_grp_id, "Time", datatype_id, time_dataspace_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
	// Clean
	H5Sclose(time_dataspace_id);
	H5Pclose(property_id);

	// get total column of the data
	obj_num = cur_output_info.size();
	total_column_num = 2; // the first two column is iteration index and time
	for (i = 0; i < obj_num; i++)
		total_column_num += cur_output_info[i].fld_num * cur_output_info[i].pcl_num;
	total_row_num = 0;
	// Get all data
	while (true)
	{
		data_len = file.getNextDataChunk(meta_handle);
		if (!data_len) break;
		data_buffer_tmp.resize(data_len);
		//std::cout << "data_buffer len: " << data_len << std::endl;
		data_buffer = (double *)data_buffer_tmp.get_buf();
		file.getData(meta_handle, data_buffer);

		chunk_row_num = data_len / sizeof(double) / total_column_num;
		total_row_num += chunk_row_num;
		
		// Output Iteration Index
		mem_dataspace_id = H5Screate_simple(1, &chunk_row_num, NULL);
		// expand time dataset in file
		H5Dset_extent(iter_index_dataset_id, &total_row_num);
		// form hyperslab
		iter_index_start = total_row_num - chunk_row_num;
		iter_index_count = chunk_row_num;
		file_dataspace_id = H5Dget_space(iter_index_dataset_id);
		H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, &iter_index_start, NULL, &iter_index_count, NULL);

		tmp_data_buf.resize(sizeof(unsigned long long) * chunk_row_num);
		tmp_data = (double *)tmp_data_buf.get_buf();
		for (data_index_i = 0; data_index_i < chunk_row_num; data_index_i++)
		{
			tmp_data[data_index_i] = data_buffer[data_index_i * total_column_num];
			//std::cout << tmp_data[i] << std::endl;
		}
		H5Dwrite(iter_index_dataset_id, ull_datatype_id, mem_dataspace_id,
					file_dataspace_id, H5P_DEFAULT, tmp_data);
		H5Sclose(mem_dataspace_id);
		H5Sclose(file_dataspace_id);

		// Output Time
		mem_dataspace_id = H5Screate_simple(1, &chunk_row_num, NULL);
		// expand time dataset in file
		H5Dset_extent(time_dataset_id, &total_row_num);
		// form hyperslab
		time_start = total_row_num - chunk_row_num;
		time_count = chunk_row_num;
		file_dataspace_id = H5Dget_space(time_dataset_id);
		H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, &time_start, NULL, &time_count, NULL);

		tmp_data_buf.resize(sizeof(double) * chunk_row_num);
		tmp_data = (double *)tmp_data_buf.get_buf();
		for (data_index_i = 0; data_index_i < chunk_row_num; data_index_i++)
		{
			tmp_data[data_index_i] = data_buffer[1 + data_index_i * total_column_num];
			//std::cout << tmp_data[i] << std::endl;
		}
		H5Dwrite(time_dataset_id, datatype_id, mem_dataspace_id, 
					file_dataspace_id, H5P_DEFAULT, tmp_data);

		H5Sclose(mem_dataspace_id);
		H5Sclose(file_dataspace_id);

		// Output data of each objects
		cur_column = 2;
		for (i = 0; i < obj_num; i++)
		{
			cur_obj_info = &cur_output_info[i];
			obj_grp_id = H5Gopen(out_grp_id, cur_obj_info->obj_name.c_str(), H5P_DEFAULT);

			fld_num = cur_obj_info->fld_num;
			pcl_num = cur_obj_info->pcl_num;

			// Output Data of each fields
			for (j = 0; j < fld_num; j++)
			{
				fld_dataset_id = H5Dopen(obj_grp_id, cur_obj_info->fld_name[j].c_str(), H5P_DEFAULT);

				dims[0] = chunk_row_num;
				dims[1] = pcl_num;
				mem_dataspace_id = H5Screate_simple(2, dims, NULL);
				// expand dataset dimension in file
				dims[0] = total_row_num;
				dims[1] = pcl_num;
				H5Dset_extent(fld_dataset_id, dims);
				// create hyperslab
				file_dataspace_id = H5Dget_space(fld_dataset_id);
				start[0] = total_row_num - chunk_row_num;
				start[1] = 0;
				count[0] = chunk_row_num;
				count[1] = pcl_num;
				H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
				
				// output data
				tmp_data_buf.resize(chunk_row_num * pcl_num * sizeof(double));
				tmp_data = (double *)tmp_data_buf.get_buf();
				for (data_index_i = 0; data_index_i < chunk_row_num; data_index_i++)
				{
					for (data_index_j = 0; data_index_j < pcl_num; data_index_j++)
					{
						tmp_data[pcl_num * data_index_i + data_index_j] = 
							data_buffer[cur_column + data_index_i * total_column_num + data_index_j];
						//std::cout << tmp_data[pcl_num * data_index_i + data_index_j] << std::endl;
					}
				}
				H5Dwrite(fld_dataset_id, datatype_id, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, tmp_data);
				cur_column += pcl_num;

				H5Sclose(mem_dataspace_id);
				H5Sclose(file_dataspace_id);

				H5Dclose(fld_dataset_id);
			}
			H5Gclose(obj_grp_id);
		}
	}

	// Add time number
	add_attr(out_grp_id, "TimeNumber", total_row_num);

	H5Tclose(ull_datatype_id);
	H5Tclose(datatype_id);
	H5Dclose(iter_index_dataset_id);
	H5Dclose(time_dataset_id);
	H5Gclose(out_grp_id);

	return 0;
}


int TmpDataToHdf5::parse_background_mesh_attribute(HTmpData meta_handle)
{
	size_t i;
	Token data_buffer_tmp;
	size_t data_len;
	char *data_buffer;
	Token tk_tmp;

	std::string key_name;

	hid_t mesh_grp_id;
	hid_t coord_dataset_id;
	hid_t coord_datatype_id;
	hid_t coord_dataspace_id;
	hsize_t coord_num;

	unsigned long long mesh_type;
	unsigned long long x_coord_num, y_coord_num;
	double x_coord, y_coord;
	std::vector<double> coord_array;

	// Assume all attribute have been written in only one chunk
	// get data
	data_len = file.getNextDataChunk(meta_handle);
	if (!data_len) return -1;
	data_buffer_tmp.resize(data_len);
	//std::cout << "data_buffer len: " << data_len << std::endl;
	data_buffer = (char *)data_buffer_tmp.get_buf();
	file.getData(meta_handle, data_buffer);

	tk_tmp.resize(128);

	mesh_grp_id = H5Gcreate(hdf5_file_id, "BackgroundMesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// mesh type
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	key_name = (char *)tk_tmp.get_buf();
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	mesh_type = *(unsigned long long *)(tk_tmp.get_buf());
	//std::cout << mesh_type << std::endl;
	add_attr(mesh_grp_id, key_name.c_str(), mesh_type);

	// mesh type name
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	add_attr(mesh_grp_id, key_name.c_str(), (char *)tk_tmp.get_buf());

	// x coord num
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name = (char *)tk_tmp.get_buf();
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	x_coord_num = *(unsigned long long *)(tk_tmp.get_buf());
	//std::cout << x_coord_num << std::endl;
	add_attr(mesh_grp_id, key_name.c_str(), x_coord_num);

	// x coord
	data_buffer = (char *)getToken(data_buffer, tk_tmp);
	//std::cout << (char *)tk_tmp.get_buf() << std::endl;
	key_name = (char *)tk_tmp.get_buf();
	for (i = 0; i < x_coord_num; i++)
	{
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		x_coord = *(double *)(tk_tmp.get_buf());
		//std::cout << x_coord << " ";
		
		coord_array.push_back(x_coord);
	}
	//std::cout << std::endl;

	// 1. Create data type (double)
	coord_datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
	H5Tset_order(coord_datatype_id, H5T_ORDER_LE);
	// 2. Create dataspace
	coord_num = x_coord_num;
	coord_dataspace_id = H5Screate_simple(1, &coord_num, NULL);
	// 3. Create dataset and write data
	coord_dataset_id = H5Dcreate(mesh_grp_id, key_name.c_str(),
		coord_datatype_id, coord_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(coord_dataset_id, coord_datatype_id, coord_dataspace_id, 
		H5S_ALL, H5P_DEFAULT, &coord_array[0]);
	coord_array.clear();
	H5Dclose(coord_dataset_id);
	// 4. Clear
	H5Tclose(coord_datatype_id);
	H5Sclose(coord_dataspace_id);

	// if mesh is two dimension
	if ((MeshType)mesh_type == MeshType::Mesh_R2D4)
	{
		// y coord num
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		y_coord_num = *(unsigned long long *)(tk_tmp.get_buf());
		//std::cout << x_coord_num << std::endl;
		add_attr(mesh_grp_id, key_name.c_str(), y_coord_num);

		// y coord
		data_buffer = (char *)getToken(data_buffer, tk_tmp);
		//std::cout << (char *)tk_tmp.get_buf() << std::endl;
		key_name = (char *)tk_tmp.get_buf();
		for (i = 0; i < y_coord_num; i++)
		{
			data_buffer = (char *)getToken(data_buffer, tk_tmp);
			y_coord = *(double *)(tk_tmp.get_buf());
			//std::cout << x_coord << " ";

			coord_array.push_back(y_coord);
		}
		//std::cout << std::endl;

		// 1. Create data type (double)
		coord_datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
		H5Tset_order(coord_datatype_id, H5T_ORDER_LE);
		// 2. Create dataspace
		coord_num = y_coord_num;
		coord_dataspace_id = H5Screate_simple(1, &coord_num, NULL);
		// 3. Create dataset and write data
		coord_dataset_id = H5Dcreate(mesh_grp_id, key_name.c_str(),
			coord_datatype_id, coord_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(coord_dataset_id, coord_datatype_id, coord_dataspace_id,
			H5S_ALL, H5P_DEFAULT, &coord_array[0]);
		coord_array.clear();
		H5Dclose(coord_dataset_id);
		// 4. Clear
		H5Tclose(coord_datatype_id);
		H5Sclose(coord_dataspace_id);
	}

	H5Gclose(mesh_grp_id);

	return 0;
}

/*

int OutputRequest::generateHDF5()
{
const char *tmp_data_file_name;
size_t tmp_data_file_name_len;
char hdf5_file_name[512];

hid_t h5_file_id;
hid_t grp_id;
hid_t dataset_id;
hid_t datatype_id;
hid_t mem_dataspace_id, file_dataspace_id;
hid_t property_id;
hid_t attribute_id;
hid_t attribute_dataspace_id;
size_t attribute_len_tmp; // used to init attributes

hsize_t dims[2];
hsize_t max_dims[2];
hsize_t chunk_dims[2];
hsize_t chunk_size;

hsize_t start[2];
hsize_t count[2];

std::string output_name;
std::string object_name;
std::string field_name;
std::string last_output_name;
std::set<std::string> output_name_set;

HTmpData meta_handle;
size_t meta_info_len;
char *meta_info_buffer;
OutputData::OutputResultInfoHeader *meta_info_header;
unsigned long long *particle_index_array;

size_t total_row_num, chunk_row_num;
size_t data_len, last_data_len;
double *data_buffer;

// Generate hdf5 file name
tmp_data_file_name = file.getFileName();
tmp_data_file_name_len = strlen(tmp_data_file_name);
memcpy(hdf5_file_name, tmp_data_file_name, tmp_data_file_name_len);
strcpy(hdf5_file_name + tmp_data_file_name_len - 8, ".h5");
// Create hdf5 file
h5_file_id = H5Fcreate(hdf5_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
if (h5_file_id < 0) return -1;

// reopen the original tmp data file
if (file.initFile(tmp_data_file_name, 0)) return -1;

// Initialize loop
grp_id = -1;
// set type of the field data
datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
H5Tset_order(datatype_id, H5T_ORDER_LE);
while (true) // for each field (metadata)
{
meta_handle = file.getNextMetaData();
if (!meta_handle) break;

// get object and field name;
parse_data_type(file.getDataType(meta_handle),output_name, object_name, field_name);

// read meta data info
meta_info_len = file.getMetaDataInfoLen(meta_handle);
meta_info_buffer = new char[meta_info_len];
file.getMetaDataInfo(meta_handle, meta_info_buffer);
meta_info_header = reinterpret_cast<OutputData::OutputResultInfoHeader *>(meta_info_buffer);

// Create or open group
if (last_output_name != output_name)
{
last_output_name = output_name;
if (output_name_set.count(output_name))
{
// Group has already been created.
if (grp_id >= 0) H5Gclose(grp_id);
grp_id = H5Dopen(h5_file_id, output_name.c_str(), H5P_DEFAULT);
}
else // if group has not been created.
{
// Create Group
output_name_set.insert(output_name);
if (grp_id >= 0) H5Gclose(grp_id);
grp_id = H5Gcreate(h5_file_id, output_name.c_str(),
H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

// ------------------- Write group attribute --------------------
// simulation type
attribute_dataspace_id = H5Screate(H5S_SCALAR);
attribute_id = H5Acreate(grp_id, "SimulationType", H5T_NATIVE_INT,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_INT, &meta_info_header->stype);
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
// object index
attribute_dataspace_id = H5Screate(H5S_SCALAR);
attribute_id = H5Acreate(grp_id, "ObjectIndex", H5T_NATIVE_ULLONG,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &meta_info_header->object_index);
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
// object name
attribute_len_tmp = strlen(object_name.c_str());
if (attribute_len_tmp)
{
attribute_dataspace_id = H5Screate_simple(1, &attribute_len_tmp, NULL);
attribute_id = H5Acreate(grp_id, "ObjectName", H5T_NATIVE_CHAR,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_CHAR, object_name.c_str());
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
}
// particle num
attribute_dataspace_id = H5Screate(H5S_SCALAR);
attribute_id = H5Acreate(grp_id, "ParticleNum", H5T_NATIVE_ULLONG,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &meta_info_header->particle_num);
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
// particle index list
particle_index_array = reinterpret_cast<unsigned long long *>
(meta_info_buffer + sizeof(OutputData::OutputResultInfoHeader));
attribute_len_tmp = meta_info_header->particle_num; // particle num
attribute_dataspace_id = H5Screate_simple(1, &attribute_len_tmp, NULL);
attribute_id = H5Acreate(grp_id, "ParticleIndexList", H5T_NATIVE_ULLONG,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_ULLONG, particle_index_array);
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
}
}

// create file data space
dims[0] = 0; // no data initially
dims[1] = meta_info_header->particle_num + 1;
max_dims[0] = H5S_UNLIMITED;
max_dims[1] = dims[1];
file_dataspace_id = H5Screate_simple(2, dims, max_dims);

// Set data set proporty
property_id = H5Pcreate(H5P_DATASET_CREATE);
H5Pset_layout(property_id, H5D_CHUNKED);
// Assume each chunk is around 50kB
chunk_size = 50 * 1024 / sizeof(double) / (meta_info_header->particle_num + 1);
chunk_size = chunk_size ? chunk_size : 1; // chunk size > 1
chunk_dims[0] = chunk_size;
chunk_dims[1] = dims[1];
H5Pset_chunk(property_id, 2, chunk_dims);

dataset_id = H5Dcreate(grp_id, field_name.c_str(), datatype_id,
file_dataspace_id, H5P_DEFAULT, property_id, H5P_DEFAULT);

H5Sclose(file_dataspace_id);
H5Pclose(property_id);

// Set properties of the dataset
// Field Type Index
attribute_dataspace_id = H5Screate(H5S_SCALAR);
attribute_id = H5Acreate(dataset_id, "FieldType", H5T_NATIVE_INT,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_INT, &meta_info_header->field);
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);
// field name
attribute_len_tmp = strlen(field_name.c_str());
attribute_dataspace_id = H5Screate_simple(1, &attribute_len_tmp, NULL);
attribute_id = H5Acreate(dataset_id, "FieldName", H5T_NATIVE_CHAR,
attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
H5Awrite(attribute_id, H5T_NATIVE_CHAR, field_name.c_str());
H5Aclose(attribute_id);
H5Sclose(attribute_dataspace_id);

// -------------------- read all data into HDF5 file -------------------
total_row_num = 0;
file.resetDataChunk(meta_handle);
data_buffer = nullptr;
last_data_len = 0;
while (true)
{
data_len = file.getNextDataChunk(meta_handle) / sizeof(double);
if (!data_len) break;

if (data_len > last_data_len)
{
if (data_buffer) delete[] data_buffer;
data_buffer = new double[data_len];
last_data_len = data_len;
}

file.getData(meta_handle, data_buffer);

chunk_row_num = data_len / (meta_info_header->particle_num + 1);
dims[0] = chunk_row_num;
dims[1] = meta_info_header->particle_num + 1;
mem_dataspace_id = H5Screate_simple(2, dims, NULL);

start[0] = total_row_num;
start[1] = 0;
count[0] = chunk_row_num;
count[1] = dims[1];

total_row_num += chunk_row_num;
dims[0] = total_row_num;
H5Dset_extent(dataset_id, dims);
file_dataspace_id = H5Dget_space(dataset_id);

H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
H5Dwrite(dataset_id, datatype_id, mem_dataspace_id,
file_dataspace_id, H5P_DEFAULT, data_buffer);

H5Sclose(mem_dataspace_id);
H5Sclose(file_dataspace_id);
}

if (data_buffer) delete[] data_buffer;
delete[] meta_info_buffer;
}

H5Gclose(grp_id);
H5Tclose(datatype_id);
H5Fclose(h5_file_id);

return 0;
}
*/