#include <iostream>
#include <cassert>

#include "ObjectOutput_1D_Hydromechanics.h"

ObjectOutput_1D_Hydromechanics::ObjectOutput_1D_Hydromechanics(
	ObjectByParticle_1D_Hydromechanics *obj) :
	ObjectOutput(obj), object(obj),
	output_funcs(nullptr), particles_ptr(nullptr), data_buffer(nullptr) {}

int ObjectOutput_1D_Hydromechanics::init(void)
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
		particles_ptr = new Particle_1D_Hydromechanics*[particle_index_num];
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

double *ObjectOutput_1D_Hydromechanics::output(double *buf)
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

bool ObjectOutput_1D_Hydromechanics::validateOutputById(unsigned long long fld_id)
{
	return getOutputFunc(fld_id) ? true : false;
}

const char *ObjectOutput_1D_Hydromechanics::outputFieldName_1D_Hydromechanics[] = {
	"x",           //0
	"n",           //1
	"stress11",    //2
	"isInMesh",    //3
	nullptr,       //4
	"mass_s",      //5
	"density_s",   //6
	"momentum1_s", //7
	"velocity1_s",  //8
	"estress11",   //9
	"strain11",    //10
	"estrain11",   //11
	"pstrain11",   //12
	nullptr,       //13
	nullptr,       //14
	"mass_f",      //15
	"density_f",   //16
	"momentum1_f", //17
	"velocity_f",  //18
	"p",           //19
};

const ObjectOutput_1D_Hydromechanics::OutputFunc 
	ObjectOutput_1D_Hydromechanics::output_fun_list[] = {
	&output_x,           //0
	&output_n,           //1
	&output_stress11,    //2
	&output_isInMesh,    //3
	nullptr,             //4
	&output_mass_s,      //5
	&output_density_s,   //6
	&output_momentum1_s, //7
	&output_velocity1_s, //8
	&output_estress11,   //9
	&output_strain11,    //10
	&output_estrain11,   //11
	&output_pstrain11,   //12
	nullptr,             //13
	nullptr,             //14
	&output_mass_f,      //15
	&output_density_f,   //16
	&output_momentum1_f, //17
	&output_velocity1_f,  //18
	&output_p,           //19
};

const size_t ObjectOutput_1D_Hydromechanics::outputFieldNum
	= sizeof(ObjectOutput_1D_Hydromechanics::output_fun_list) /
	  sizeof(ObjectOutput_1D_Hydromechanics::output_fun_list[0]);

ObjectOutput_1D_Hydromechanics::OutputFunc
ObjectOutput_1D_Hydromechanics::getOutputFunc(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? output_fun_list[fld_id] : nullptr;
}

const char *ObjectOutput_1D_Hydromechanics::
	getFieldName(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? 
		outputFieldName_1D_Hydromechanics[fld_id] : nullptr;
}

void ObjectOutput_1D_Hydromechanics::output_x()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->x;
}

void ObjectOutput_1D_Hydromechanics::output_n()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->n;
}

void ObjectOutput_1D_Hydromechanics::output_stress11()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress11;
}

void ObjectOutput_1D_Hydromechanics::output_isInMesh()
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

void ObjectOutput_1D_Hydromechanics::output_mass_s()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->mass_s;
}

void ObjectOutput_1D_Hydromechanics::output_density_s()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->density_s;
}

void ObjectOutput_1D_Hydromechanics::output_momentum1_s()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1_s;
}

void ObjectOutput_1D_Hydromechanics::output_velocity1_s()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1_s / particles_ptr[i]->mass_s;
}

void ObjectOutput_1D_Hydromechanics::output_estress11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estress11;
}

void ObjectOutput_1D_Hydromechanics::output_strain11()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->strain11;
}

void ObjectOutput_1D_Hydromechanics::output_estrain11()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estrain11;
}

void ObjectOutput_1D_Hydromechanics::output_pstrain11()
{
	size_t i;
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->pstrain11;
}

void ObjectOutput_1D_Hydromechanics::output_mass_f()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->mass_f;
}

void ObjectOutput_1D_Hydromechanics::output_density_f()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->density_f;
}

void ObjectOutput_1D_Hydromechanics::output_momentum1_f()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1_f;
}

void ObjectOutput_1D_Hydromechanics::output_velocity1_f()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1_f / particles_ptr[i]->mass_f;
}

void ObjectOutput_1D_Hydromechanics::output_p()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->p;
}
