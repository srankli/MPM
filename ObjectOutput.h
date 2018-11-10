#ifndef __OBJECTOUTPUT_H__
#define __OBJECTOUTPUT_H__

#include <string>
#include <vector>

#include "Particle.h"

class Output;

// Each Output class contains multiple ObjectOutput class
class ObjectOutput
{
	friend Output;
protected:
	ObjectByParticle *obj_base;
	// Field type
	size_t field_type_num;
	std::vector<unsigned long long> field_type_array;
	// Particle index
	size_t particle_index_num;
	std::vector<size_t> particle_index_array;
	// if all_particles is true, means include all particles of an object
	bool all_particles;
	// number of output
	size_t output_num;
public:
	/* If pcls == nullptr, it means use all particles of the object. */
	ObjectOutput(ObjectByParticle *obj) :
		obj_base(obj), output_num(0),
		field_type_num(0), particle_index_num(0), all_particles(false) {}
	virtual ~ObjectOutput() = 0
	{
		field_type_array.clear();
		particle_index_array.clear();
	}
	int addFields(std::vector<unsigned long long> *flds)
	{
		size_t i, num_tmp;

		// count valid field index
		num_tmp = flds->size();
		for (i = 0; i < num_tmp; i++)
		{
			if (validateOutputById((*flds)[i]))
				field_type_array.push_back((*flds)[i]);
		}
		field_type_num = field_type_array.size();
		return 0;
	}

	int addParticles(std::vector<size_t> *pcls = nullptr)
	{
		size_t i, num_tmp;
		Particle *pcl_iter;

		if (!pcls) // all paricles of the object has been selected
		{
			all_particles = true;
			particle_index_num = obj_base->getParticleNum();
			pcl_iter = obj_base->getFirstParticle();
			for (i = 0; i < particle_index_num; i++)
			{
				particle_index_array.push_back(pcl_iter->index);
				//std::cout << pcl_iter->index << std::endl;
				pcl_iter = obj_base->getNextParticle(pcl_iter);
			}
		}
		else
		{
			// count valid particle index
			num_tmp = pcls->size();
			for (i = 0; i < num_tmp; i++)
			{
				if (obj_base->validateParticleId((*pcls)[i]))
					particle_index_array.push_back((*pcls)[i]);
			}
			particle_index_num = particle_index_array.size();
		}
		return 0;
	}


	inline SimulationType getSimulationType(void) noexcept { return obj_base->getSimulationType(); }
	inline size_t getObjectIndex(void) noexcept { return obj_base->getIndex(); }
	inline const std::string &getObjectName() noexcept { return obj_base->getName(); }
	inline size_t getOutputNum(void) noexcept { return output_num; }

	// Initialzie this output
	virtual int init(void) = 0;
	// Output data to buffer
	// return current location of the buffer
	virtual double *output(double *buffer) = 0;
	// check if the fld_id has corresponding output
	virtual bool validateOutputById(unsigned long long fld_id) = 0;
	virtual const char *getFieldName(size_t fld_id) = 0;
};

#endif