#ifndef __OBJECTOUTPUT_1D_HYDROMECHANICS_H__
#define __OBJECTOUTPUT_1D_HYDROMECHANICS_H__

#include "ObjectOutput.h"
#include "Particle_1D_Hydromechanics.h"

class OutputRequest;

enum class OutputFieldType_1D_Hydromechanics : unsigned long long
{
	x = 0,
	n = 1,
	stress11 = 2,
	isInMesh = 3,

	mass_s = 5,
	density_s = 6,
	momentum1_m = 7,
	velocity1_s = 8,
	estress11 = 9,
	strain11 = 10,
	estrain11 = 11,
	pstrain11 = 12,

	mass_f = 15,
	density_f = 16,

	velocity1_f = 18,
	p = 19
};

class ObjectOutput_1D_Hydromechanics : public ObjectOutput
{
protected:
	// output function table
	typedef void(ObjectOutput_1D_Hydromechanics::*OutputFunc)(void);

	ObjectByParticle_1D_Hydromechanics *object;
	// Output function list
	OutputFunc *output_funcs;
	// Particle pointers list
	Particle_1D_Hydromechanics **particles_ptr;

	// data buffer, needed to be updated by every output(double*)
	double *data_buffer;

public:
	ObjectOutput_1D_Hydromechanics(ObjectByParticle_1D_Hydromechanics *obj);
	~ObjectOutput_1D_Hydromechanics()
	{
		if (output_funcs) delete[] output_funcs;
		if (particles_ptr) delete[] particles_ptr;
	}

	int init(void);
	double *output(double *buf);
	bool validateOutputById(unsigned long long fld_id);

protected:
	// output function table
	static const size_t outputFieldNum;
	static const OutputFunc output_fun_list[];
	OutputFunc getOutputFunc(unsigned long long fld_id);

	void output_x();
	void output_n();
	void output_stress11();
	void output_isInMesh();
	void output_mass_s();
	void output_density_s();
	void output_momentum1_m();
	void output_velocity1_s();
	void output_estress11();
	void output_strain11();
	void output_estrain11();
	void output_pstrain11();
	void output_mass_f();
	void output_density_f();
	void output_velocity1_f();
	void output_p();

	// field name string
	static const char *outputFieldName_1D_Hydromechanics[];
public:
	const char *getFieldName(unsigned long long fld_id);
};

#endif