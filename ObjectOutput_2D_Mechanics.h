#ifndef __OBJECTOUTPUT_2D_MECHANICS_H__
#define __OBJECTOUTPUT_2D_MECHANICS_H__

#include "ObjectOutput.h"
#include "Particle_2D_Mechanics.h"

class OutputRequest;

enum class OutputFieldType_2D_Mechanics : unsigned long long
{
	x = 0,
	y = 1,
	isInMesh = 3,

	mass = 5,
	density = 6,
	momentum1 = 7,
	momentum2 = 8,
	velocity1 = 9,
	velocity2 = 10,

	stress11 = 11,
	stress22 = 12,
	stress33 = 13,
	stress12 = 14,
	stress23 = 15,
	stress31 = 16,

	strain11 = 20,
	strain22 = 21,
	strain12 = 22,
	estrain11 = 23,
	estrain22 = 24,
	estrain12 = 25,
	pstrain11 = 26,
	pstrain22 = 27,
	pstrain12 = 28
};

class ObjectOutput_2D_Mechanics : public ObjectOutput
{
protected:
	// output function table
	typedef void(ObjectOutput_2D_Mechanics::*OutputFunc)(void);

	ObjectByParticle_2D_Mechanics *object;
	// Output function list
	OutputFunc *output_funcs;
	// Particle pointers list
	Particle_2D_Mechanics **particles_ptr;

	// data buffer, needed to be updated by every output(double*)
	double *data_buffer;

public:
	ObjectOutput_2D_Mechanics(ObjectByParticle_2D_Mechanics *obj);
	~ObjectOutput_2D_Mechanics()
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
	void output_y();
	void output_isInMesh();
	void output_mass();
	void output_density();
	void output_momentum1();
	void output_momentum2();
	void output_velocity1();
	void output_velocity2();
	void output_stress11();
	void output_stress22();
	void output_stress33();
	void output_stress12();
	void output_stress23();
	void output_stress31();
	void output_strain11();
	void output_strain22();
	void output_strain12();
	void output_estrain11();
	void output_estrain22();
	void output_estrain12();
	void output_pstrain11();
	void output_pstrain22();
	void output_pstrain12();

	// field name string
	static const char *outputFieldName_2D_Mechanics[];
public:
	const char *getFieldName(unsigned long long fld_id);
};

#endif