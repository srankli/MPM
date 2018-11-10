#ifndef __TMPDATATOHDF5_H__
#define __TMPDATATOHDF5_H__

#include <string>
#include <map>

#define H5_BUILT_AS_DYNAMIC_LIB
#include "hdf5.h"
#pragma comment(lib, "hdf5.lib")

#include "FileBuffer.h"

class TmpDataToHdf5
{
protected:
	std::string tmp_data_file_name;
	std::string hdf5_file_name;
	TmpDataFile file;

public:
	TmpDataToHdf5();
	~TmpDataToHdf5();

	int init(const char *file_name);
	int init(OutputRequest &out_req);
	int generate(void);
	void complete(void);

	// assume that all attribute has been written in one chunk
	// meta data of output data must directly follow attribute in TmpDataFile 
	int parse_output_attribute(HTmpData meta_handle);
	int parse_output_data(HTmpData meta_handle);

	int parse_background_mesh_attribute(HTmpData meta_handle);

protected:
	//std::map<unsigned long long, hid_t> output_id_map;
	hid_t hdf5_file_id;
	struct ObjectInfo
	{
		size_t obj_index;
		std::string obj_name;
		size_t fld_num;
		std::vector<unsigned long long> fld_type;
		std::vector<std::string> fld_name;
		size_t pcl_num;
	};
	// output_info
	hid_t cur_output_id;
	std::string cur_output_name;
	std::vector<ObjectInfo> cur_output_info;
};

#endif