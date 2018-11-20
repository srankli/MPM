#include "Particle_1D_Hydromechanics.h"

int ObjectByParticle_1D_Hydromechanics::addParticle(
	ParticleParam_1D_Hydromechanics *pcl_param,
	ConstitutiveModelParam *pcl_cm)
{
	Particle_1D_Hydromechanics *pcl;
	double volume_s, volume_f;

	if (!pcl_param || !pcl_cm)
		return -1;

	++particleNum;

	pcl = static_cast<Particle_1D_Hydromechanics *>(particles_mem.alloc());
	pcl->index = ++curParticleIndex;
	pcl->object = this;

	// --------- Solid Phase ---------
	pcl->x = pcl_param->x;
	pcl->mass_s = pcl_param->mass_s;
	pcl->density_s = pcl_param->density_s;
	pcl->momentum1_m = pcl_param->momentum1_m;
	pcl->estress11 = pcl_param->estress11;
	pcl->strain11 = pcl_param->strain11;
	pcl->estrain11 = pcl_param->estrain11;
	pcl->pstrain11 = pcl_param->pstrain11;

	// --------- Fluid Phase ---------
	pcl->p = pcl_param->p;
	pcl->mass_f = pcl_param->mass_f;
	pcl->density_f = pcl_param->density_f;
	pcl->Kf = pcl_param->Kf;
	// --- Solid - Fluid Interaction ---
	pcl->k = pcl_param->k; // permeability (L/T)
	
	volume_s = pcl_param->mass_s / pcl_param->density_s;
	volume_f = pcl_param->mass_f / pcl_param->density_f;
	pcl->n = volume_f / (volume_s + volume_f);
	pcl->stress11 = pcl_param->estress11 - pcl_param->p;

	constitutiveModels_mem.add_model(pcl_cm);

	return 0;
}

void ObjectByParticle_1D_Hydromechanics::finish_init(void)
{
	ObjectByParticle::finish_init();
	particles = static_cast<Particle_1D_Hydromechanics *>(particles_mem.get_mem());
}
