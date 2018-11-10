#include <iostream>
#include <cassert>

#include "ObjectOutput_2D_Mechanics.h"

ObjectOutput_2D_Mechanics::ObjectOutput_2D_Mechanics(
	ObjectByParticle_2D_Mechanics *obj) :
	ObjectOutput(obj), object(obj),
	output_funcs(nullptr), particles_ptr(nullptr), data_buffer(nullptr) {}


int ObjectOutput_2D_Mechanics::init(void)
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
		particles_ptr = new Particle_2D_Mechanics*[particle_index_num];
		if (!particles_ptr)	return -1;
		for (i = 0; i < particle_index_num; i++)
			particles_ptr[i] = object->getParticleById(particle_index_array[i]);
	}

	output_num = field_type_num * particle_index_num;

	return 0;
}


double *ObjectOutput_2D_Mechanics::output(double *buf)
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

bool ObjectOutput_2D_Mechanics::validateOutputById(unsigned long long fld_id)
{
	return getOutputFunc(fld_id) ? true : false;
}

// Field name list
const char *ObjectOutput_2D_Mechanics::outputFieldName_2D_Mechanics[] = {
	"x",         //0
	"y",         //1
	nullptr,     //2
	"isInMesh",  //3
	nullptr,     //4
	"mass",      //5
	"density",   //6
	"momentum1", //7
	"momentum2", //8
	"velocity1", //9
	"velocity2", //10
	"stress11",  //11
	"stress22",  //12
	"stress33",  //13
	"stress12",  //14
	"stress23",  //15
	"stress31",  //16
	nullptr,     //17
	nullptr,     //18
	nullptr,     //19
	"strain11",  //20
	"strain22",  //21
	"strain12",  //22
	"estrain11", //23
	"estrain22", //24
	"estrain12", //25
	"pstrain11", //26
	"pstrain22", //27
	"pstrain12", //28
};

const ObjectOutput_2D_Mechanics::OutputFunc 
	ObjectOutput_2D_Mechanics::output_fun_list[] = {
	&output_x,          //0
	&output_y,          //1
	nullptr,            //2
	&output_isInMesh,   //3
	nullptr,            //4
	&output_mass,       //5
	&output_density,    //6
	&output_momentum1,  //7
	&output_momentum2,  //8
	&output_velocity1,  //9
	&output_velocity2,  //10
	&output_stress11,   //11
	&output_stress22,   //12
	&output_stress33,   //13
	&output_stress12,   //14
	&output_stress23,   //15
	&output_stress31,   //16
	nullptr,            //17
	nullptr,            //18
	nullptr,            //19
	&output_strain11,   //20
	&output_strain22,   //21
	&output_strain12,   //22
	&output_estrain11,  //23
	&output_estrain22,  //24
	&output_estrain12,  //25
	&output_pstrain11,  //26
	&output_pstrain22,  //27
	&output_pstrain12,  //28
};

const size_t ObjectOutput_2D_Mechanics::outputFieldNum
	= sizeof(ObjectOutput_2D_Mechanics::output_fun_list) /
	  sizeof(ObjectOutput_2D_Mechanics::output_fun_list[0]);

ObjectOutput_2D_Mechanics::OutputFunc
	ObjectOutput_2D_Mechanics::getOutputFunc(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? output_fun_list[fld_id] : nullptr;
}

const char *ObjectOutput_2D_Mechanics::getFieldName(unsigned long long fld_id)
{
	return fld_id < outputFieldNum ? outputFieldName_2D_Mechanics[fld_id] : nullptr;
}


void ObjectOutput_2D_Mechanics::output_x()
{
	size_t i;
	// output each fields
	for (i = 0; i < particle_index_num; i++)
	{
		data_buffer[i] = particles_ptr[i]->x;
		//std::cout << data_buffer[i] << std::endl;
	}
	//std::cout << particles_ptr[0]->x << std::endl;
}

void ObjectOutput_2D_Mechanics::output_y()
{
	size_t i;
	// output each fields
	for (i = 0; i < particle_index_num; i++)
	{
		data_buffer[i] = particles_ptr[i]->y;
		//std::cout << data_buffer[i] << std::endl;
	}
}

void ObjectOutput_2D_Mechanics::output_isInMesh()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
	{
		if (particles_ptr[i]->isInMesh)
			data_buffer[i] = 1.0;
		else
			data_buffer[i] = 0.0;
	}
}

void ObjectOutput_2D_Mechanics::output_mass()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->mass;
}

void ObjectOutput_2D_Mechanics::output_density()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->density;
}

void ObjectOutput_2D_Mechanics::output_momentum1()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1;
}

void ObjectOutput_2D_Mechanics::output_momentum2()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum2;
}

void ObjectOutput_2D_Mechanics::output_velocity1()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum1 / particles_ptr[i]->mass;
}

void ObjectOutput_2D_Mechanics::output_velocity2()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->momentum2 / particles_ptr[i]->mass;
	//std::cout << data_buffer[5] << std::endl;
}

void ObjectOutput_2D_Mechanics::output_stress11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress11;
}

void ObjectOutput_2D_Mechanics::output_stress22()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress22;
}

void ObjectOutput_2D_Mechanics::output_stress33()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress33;
}

void ObjectOutput_2D_Mechanics::output_stress12()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress12;
}

void ObjectOutput_2D_Mechanics::output_stress23()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress23;
}

void ObjectOutput_2D_Mechanics::output_stress31()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->stress31;
}

void ObjectOutput_2D_Mechanics::output_strain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->strain11;

}

void ObjectOutput_2D_Mechanics::output_strain22()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->strain22;

}

void ObjectOutput_2D_Mechanics::output_strain12()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->strain12;

}

void ObjectOutput_2D_Mechanics::output_estrain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estrain11;

}

void ObjectOutput_2D_Mechanics::output_estrain22()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estrain22;

}

void ObjectOutput_2D_Mechanics::output_estrain12()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->estrain12;

}

void ObjectOutput_2D_Mechanics::output_pstrain11()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->pstrain11;
}

void ObjectOutput_2D_Mechanics::output_pstrain22()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->pstrain22;
}

void ObjectOutput_2D_Mechanics::output_pstrain12()
{
	size_t i;
	//  output each fields
	for (i = 0; i < particle_index_num; i++)
		data_buffer[i] = particles_ptr[i]->pstrain12;
}