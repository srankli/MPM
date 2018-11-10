#ifndef __OUTPUTREQUEST_H__
#define __OUTPUTREQUEST_H__

#include <vector>
#include <string>
#include <cstdlib>

#include "FileBuffer.h"

// ***Include all OutputData of all simulation type here
// Remeber also to modify addOutput functions when adding new type of output
#include "ObjectOutput.h"
#include "ObjectOutput_1D_Mechanics.h"
//#include "OutputData_1D_Hydromechanics.h"
#include "ObjectOutput_2D_Mechanics.h"

class OutputRequest;
// used for mesh output
class Mesh;
class TmpDataToHdf5;

class FileDataKey
{
public:
	struct KeyType
	{
		const char *key;
		size_t len;
		KeyType(const char *key_string) :
			key(key_string), len(strlen(key_string)) {}
	};
	static KeyType OutputData;
	static KeyType OutputAttribute;
	static KeyType OutputIndex;
	static KeyType OutputName;
	static KeyType ObjectNumber;
	static KeyType ObjectIndex;
	static KeyType ObjectName;
	static KeyType SimulationType;
	static KeyType SimulationTypeName;
	static KeyType FieldNumber;
	static KeyType FieldType;
	static KeyType FieldName;
	static KeyType ParticleNumber;
	static KeyType ParticleIndex;
	
	static KeyType BackgroundMeshAttribute;
	static KeyType MeshType;
	static KeyType MeshTypeName;
	static KeyType XCoordNum;
	static KeyType XCoord;
	static KeyType YCoordNum;
	static KeyType YCoord;
};


class ObjectOutputManager : public FlexibleSizeMemory
{
public:
	ObjectOutputManager() {}
	~ObjectOutputManager();

	// ***Need to modify this function when new simulation is added
	int add_object_output(ObjectByParticle *obj,
		std::vector<unsigned long long> *flds, std::vector<size_t> *pcls = nullptr);
	
	inline ObjectOutput *get_first() noexcept
	{
		return static_cast<ObjectOutput *>(FlexibleSizeMemory::get_first());
	}
	inline ObjectOutput *get_next(ObjectOutput *prev_model) noexcept
	{
		return static_cast<ObjectOutput *>(FlexibleSizeMemory::get_next(prev_model));
	}
};


class Output
{
	friend OutputRequest;
protected:
	// used to initialized index
	static size_t curIndex;
	size_t index;
	std::string name;
	Output *next;

	size_t output_num;
	double t_time;
	double t_inc;
	double t_tol;

	ObjectOutputManager output_mem;
	size_t object_num;

	AttributeBuffer attr_buffer;

	size_t buffer_size_recomended_value; // default 10k
	DataBuffer<double> data_buffer;

public:
	/*  Parameter:
	1. out_name: name of the output
	2. num: times of output
	3. bsrv: size of buffer
	*/
	Output(const char *out_name, size_t num, Output *nt, size_t bsrv);
	~Output() {}

	inline size_t getIndex() noexcept { return index; }
	inline const std::string &getName() noexcept { return name; }

	int addObject(ObjectByParticle *obj,
		std::vector<unsigned long long> *flds,
		std::vector<size_t> *pcls = nullptr)
	{
		int res = output_mem.add_object_output(obj, flds, pcls);
		if (!res) ++object_num;
		return res;
	}
	int init(TmpDataFile &outFile, double step_time);
	double output(unsigned long long iter_index, double t_cur, bool mustOutput = false);
	void complete(void);
};


class OutputRequest
{
	friend TmpDataToHdf5;
protected:
	double step_time;
	std::string file_name;

	double next_output_time; // used by output(t_cur)
	TmpDataFile file;

	Output *top;
	Output *tail;
	size_t output_num;
public:
	OutputRequest() : step_time(0.0), next_output_time(0.0),
		top(nullptr), tail(nullptr), output_num(0) {}
	~OutputRequest();

	inline void setFileName(const char *fname) noexcept { file_name = fname; }
	inline const std::string &getFileName(void) noexcept { return file_name; }
	inline void setStepTime(const double st) noexcept {	step_time = st;	}
	inline double getStepTime(void) noexcept { return step_time; }
	inline size_t getOutputNum(void) noexcept { return output_num; }
	Output *addOutput(const char *out_name,
		size_t output_num, /* times of output in this step */
		size_t buf_size = 10240 /* default 10k */);
	int finish_init(void);
	void output(unsigned long long iter_index, double t_cur, bool mustOutput = false);
	void complete(void);
	void reset(void);

	void outputMeshGeometry(Mesh *mh);

public:
	// convert current TmpDataFile *file to HDF5
	//int generateHDF5();
};

#endif