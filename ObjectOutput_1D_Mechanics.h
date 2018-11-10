#ifndef __OBJECTOUTPUT_1D_MECHANICS_H__
#define __OBJECTOUTPUT_1D_MECHANICS_H__

#include "ObjectOutput.h"
#include "Particle_1D_Mechanics.h"

class OutputRequest;

enum class OutputFieldType_1D_Mechanics : unsigned long long
{
	x = 0,
	isInMesh = 3,

	mass = 5,
	density = 6,
	momentum1 = 7,
	velocity1 = 8,
	stress11 = 9,
	strain11 = 10,
	estrain11 = 11,
	pstrain11 = 12,
};


class ObjectOutput_1D_Mechanics : public ObjectOutput
{
protected:
	// output function table
	typedef void(ObjectOutput_1D_Mechanics::*OutputFunc)(void);

	ObjectByParticle_1D_Mechanics *object;
	// Output function list
	OutputFunc *output_funcs;
	// Particle pointers list
	Particle_1D_Mechanics **particles_ptr;
	
	// data buffer, needed to be updated by every output(double*)
	double *data_buffer;

public:
	ObjectOutput_1D_Mechanics(ObjectByParticle_1D_Mechanics *obj);
	~ObjectOutput_1D_Mechanics()
	{
		if (output_funcs) delete[] output_funcs;
		if (particles_ptr) delete[] particles_ptr;
	}

	int init(void);
	double *output(double *buf);
	bool validateOutputById(unsigned long long fld_id);

protected:
	// number of fields that are available for output 
	static const size_t outputFieldNum;
	static const OutputFunc output_fun_list[];
	OutputFunc getOutputFunc(unsigned long long fld_id);

	void output_x();
	void output_isInMesh();
	void output_mass();
	void output_density();
	void output_momentum1();
	void output_velocity1();
	void output_stress11();
	void output_strain11();
	void output_estrain11();
	void output_pstrain11();

	// field name string
	static const char *outputFieldName_1D_Mechanics[];
public:
	const char *getFieldName(unsigned long long fld_id);
};

#endif