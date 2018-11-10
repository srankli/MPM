#include <iostream>
#include <cassert>

#include "ObjectOutput_1D_Mechanics.h"

ObjectOutput_1D_Mechanics::ObjectOutput_1D_Mechanics(
	ObjectByParticle_1D_Mechanics *obj) :
	ObjectOutput(obj), object(obj),
	output_funcs(nullptr), particles_ptr(nullptr), data_buffer(nullptr) {}

int ObjectOutput_1D_Mechanics::init(void)
{
	size_t i;

	// Initialize field output function list
	if (field_type_num)
	{
		output_funcs = new OutputFunc[field_type_num];
		if (!output_funcs)	return -1;
		for (i = 0; i < field_type_num; i++)
			output_funcs[i] = getOutputFunc(field_type_array[i]);
	}

	// Initialize particle pointer list
	if (particle_index_num)
	{
		particles_ptr = new Particle_1D_Mechanics*[particle_index_num];
		if (!particles_ptr)	return -1;
		for (i = 0; i < particle_index_num; i++)
		{
			particles_ptr[i] = object->getParticleById(particle_index_array[i]);
			//std::cout << particle_index_array[i] << std::endl;
		}
	}

	output_num = field_type_num * particle_index_num;

	return 0;
}

double *ObjectOutput_1D_Mechanics::output(double *buf)
{
	size_t i;

	data_buffer = buf;
	for (i = 0; i < field_type_num; i++)
	{
		(this->*output_funcs[i])();
		data_buffer += particle_index_num;
	}
	return data_buffer;
}

bool ObjectOutput_1D_Mechanics::validateOutputById(unsigned long long fld_id)
{
	return getOutputFunc(fld_id) ? true : false;
}


// Field name list
const char *ObjectOutput_1D_Mechanics::outputFieldName_1D_Mechanics[] = {
	"x",           //0
	nullptr,       //1
	nullptr,       //2
	"isInMesh",    //3
	nullptr,       //4
	"mass",      //5
	"density",   //6
	"momentum1", //7
	"velocity1", //8
	"stress11",    //9
	"strain11",    //10
	"estrain11",   //11
	"pstrain11",   //12
};

// Function pointers list
const ObjectOutput_1D_Mechanics::OutputFunc
	ObjectOutput_1D_Mechanics::output_fun_list[] = {
	&output_x,          //0
	nullptr,            //1
	nullptr,            //2
	&output_isInMesh,   //3
	nullptr,            //4
	&output_mass,       //5
	&output_density,    //6
	&output_momentum1,  //7
	&output_velocity1,  //8
	&output_stress11,   //9
	&output_strain11,   //10
	&output_estrain11,  //11
	&output_pstrain11,  //12
};

const size_t ObjectOutput_1D_Mechanics::outputFieldNum
	= sizeof(ObjectOutput_1D_Mechanics::output_fun_list) /
	  sizeof(ObjectOutput_1D_Mechanics::output_fun_list[0]);

ObjectOutput_1D_Mechanics::OutputFunc 
	ObjectOutput_1D_Mechanics::getOutputFunc(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? output_fun_list[fld_id] : nullptr;
}

const char *ObjectOutput_1D_Mechanics::getFieldName(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? outputFieldName_1D_Mechanics[fld_id] : nullptr;
}


void ObjectOutput_1D_Mechanics::output_x()
{
	size_t i;
	// output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->x;
}

void ObjectOutput_1D_Mechanics::output_isInMesh()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
	{
		if (particles_ptr[i]->isInMesh)
			data_buffer[i] = 1.0;
		else
			data_buffer[i] = 0.0;
	}
}

void ObjectOutput_1D_Mechanics::output_mass()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->mass;
}

void ObjectOutput_1D_Mechanics::output_density()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->density;
}

void ObjectOutput_1D_Mechanics::output_momentum1()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1;
}

void ObjectOutput_1D_Mechanics::output_velocity1()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1 / particles_ptr[i]->mass;
}

void ObjectOutput_1D_Mechanics::output_stress11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress11;
}

void ObjectOutput_1D_Mechanics::output_strain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->strain11;

}

void ObjectOutput_1D_Mechanics::output_estrain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estrain11;

}

void ObjectOutput_1D_Mechanics::output_pstrain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->pstrain11;
}
