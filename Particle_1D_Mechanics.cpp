#include "Particle_1D_Mechanics.h"

int ObjectByParticle_1D_Mechanics::addParticle(
	ParticleParam_1D_Mechanics *pcl_param, ConstitutiveModelParam *pcl_cm)
{
	Particle_1D_Mechanics *pcl;

	if (!pcl_param || !pcl_cm)
		return -1;

	++particleNum;

	pcl = static_cast<Particle_1D_Mechanics *>(particles_mem.alloc());
	pcl->index = ++curParticleIndex;
	pcl->object = this;
	pcl->x = pcl_param->x;
	pcl->mass = pcl_param->mass;
	pcl->density = pcl_param->density;
	pcl->momentum1 = pcl_param->momentum1;
	pcl->stress11 = pcl_param->stress11;
	pcl->strain11 = pcl_param->strain11;
	pcl->estrain11 = pcl_param->estrain11;
	pcl->pstrain11 = pcl_param->pstrain11;
	
	constitutiveModels_mem.add_model(pcl_cm);

	return 0;
}


void ObjectByParticle_1D_Mechanics::finish_init()
{
	ObjectByParticle::finish_init();
	particles = static_cast<Particle_1D_Mechanics *>(particles_mem.get_mem());

	/*
	// print the whole object
	size_t i;
	std::cout << "Particle info:" << std::endl;
	for (i = 0; i < particleNum; i++)
	{
		std::cout << "No.: " << particles[i].index << " CM: "
			<< (unsigned int)(particles[i].cm->getType()) << std::endl;
	}
	std::cout << "Mass force BC:" << std::endl;
	for (i = 0; i < massForceBCNum; i++)
	{
		std::cout << "Pcl_id:" << massForceBCs[i].particle->index << " Value: "
			<< massForceBCs[i].bodyForce1 << std::endl;
	}
	std::cout << "Surface force BC:" << std::endl;
	for (i = 0; i < surfaceForceBCNum; i++)
	{
		std::cout << "Pcl_id:" << surfaceForceBCs[i].particle->index << " Value: "
			<< surfaceForceBCs[i].surfaceForce1 << std::endl;
	}
	*/
}