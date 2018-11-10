#include "Particle_1D_Hydromechanics.h"

int ObjectByParticle_1D_Hydromechanics::addParticle(
	Particle_1D_Hydromechanics *pcl_param,
	ConstitutiveModelParam *pcl_cm)
{
	Particle_1D_Hydromechanics *pcl;

	if (!pcl_param || !pcl_cm)
		return -1;

	++particleNum;

	pcl = static_cast<Particle_1D_Hydromechanics *>(particles_mem.alloc());
	pcl->index = ++curParticleIndex;
	pcl->object = this;
	// --------- Solid Phase ---------
	pcl->x = pcl_param->x;
	pcl->n = pcl_param->n;
	pcl->stress11 = pcl_param->stress11;
	pcl->mass_s = pcl_param->mass_s;
	pcl->density_s = pcl_param->density_s;
	pcl->momentum1_s = pcl_param->momentum1_s;
	pcl->estrain11 = pcl_param->estrain11;
	pcl->strain11 = pcl_param->strain11;
	pcl->estrain11 = pcl_param->estrain11;
	pcl->pstrain11 = pcl_param->pstrain11;
	// --------- Fluid Phase ---------
	pcl->p = pcl_param->p;
	pcl->mass_f = pcl_param->mass_f;
	pcl->density_f = pcl_param->density_f;
	pcl->momentum1_f = pcl_param->momentum1_f;
	pcl->Kf = pcl_param->Kf;
	// --- Solid - Fluid Interaction ---
	pcl->kappa = pcl_param->kappa; // permeability
	
	constitutiveModels_mem.add_model(pcl_cm);

	return 0;
}

void ObjectByParticle_1D_Hydromechanics::finish_init(void)
{
	ObjectByParticle::finish_init();
	particles = static_cast<Particle_1D_Hydromechanics *>(particles_mem.get_mem());
}
