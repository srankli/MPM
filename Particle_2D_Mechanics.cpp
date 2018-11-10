#include "Particle_2D_Mechanics.h"

int ObjectByParticle_2D_Mechanics::addParticle(
	ParticleParam_2D_Mechanics *pcl_param, ConstitutiveModelParam *pcl_cm)
{
	Particle_2D_Mechanics *pcl;

	if (!pcl_param || !pcl_cm)
		return -1;

	++particleNum;

	pcl = new(particles_mem.alloc())Particle_2D_Mechanics;
	pcl->index = ++curParticleIndex;
	pcl->object = this;

	pcl->x = pcl_param->x;
	pcl->y = pcl_param->y;
	pcl->mass = pcl_param->mass;
	pcl->density = pcl_param->density;

	pcl->momentum1 = pcl_param->momentum1;
	pcl->momentum2 = pcl_param->momentum2;

	pcl->stress11 = pcl_param->stress11;
	pcl->stress22 = pcl_param->stress22;
	pcl->stress33 = pcl_param->stress33;
	pcl->stress12 = pcl_param->stress12;
	pcl->stress23 = pcl_param->stress23;
	pcl->stress31 = pcl_param->stress31;

	pcl->strain11 = pcl_param->strain11;
	pcl->strain22 = pcl_param->strain22;
	pcl->strain12 = pcl_param->strain12;

	pcl->estrain11 = pcl_param->estrain11;
	pcl->estrain22 = pcl_param->estrain22;
	pcl->estrain12 = pcl_param->estrain12;

	pcl->pstrain11 = pcl_param->pstrain11;
	pcl->pstrain22 = pcl_param->pstrain22;
	pcl->pstrain12 = pcl_param->pstrain12;

	constitutiveModels_mem.add_model(pcl_cm);

	return 0;
}

void ObjectByParticle_2D_Mechanics::finish_init()
{
	ObjectByParticle::finish_init();
	particles = static_cast<Particle_2D_Mechanics *>(particles_mem.get_mem());

	/*
	// print the whole object
	size_t i;
	std::cout << "Particle info:" << std::endl;
	for (i = 0; i < particleNum; i++)
	{
		std::cout << "No.: " << particles[i].index << " CM: "
			<< (unsigned int)(particles[i].cm->getType()) << std::endl;
	}
	if (massForceBCNum) std::cout << "Mass force BC:" << std::endl;
	for (i = 0; i < massForceBCNum; i++)
	{
		std::cout << "Pcl_id:" << massForceBCs[i].particle->index << " Value: "
			<< massForceBCs[i].bodyForce1 << std::endl;
	}
	if (surfaceForceBCNum) std::cout << "Surface force BC:" << std::endl;
	for (i = 0; i < surfaceForceBCNum; i++)
	{
		std::cout << "Pcl_id:" << surfaceForceBCs[i].particle->index << " Value: "
			<< surfaceForceBCs[i].surfaceForce1 << std::endl;
	}
	*/
}

void ObjectByParticle_2D_Mechanics::addParticleFromTriElement(
	double x1, double y1,
	double x2, double y2,
	double x3, double y3,
	ParticleParam_2D_Mechanics *pcl_param,
	ConstitutiveModelParam *pcl_cm,
	double max_characteristic_size)
{
	size_t i, j;
	size_t div_num;
	double element_area;
	double L1, L2, L3;
	double sub_x1, sub_y1;
	double sub_x2, sub_y2;
	double sub_x3, sub_y3;
	double sub_elem_area;

	element_area = abs(x1 * y2 + x2 * y3 + x3 * y1
		- x1 * y3 - x2 * y1 - x3 * y2) / 2.0;
	if (max_characteristic_size < 0)
		div_num = 1;
	else
		div_num = int(sqrt(element_area / 3.0) / max_characteristic_size) + 1;
	sub_elem_area = element_area / (div_num * div_num);

	for (i = 0; i < div_num; i++)
	{
		for (j = 0; j <= i; j++)
		{
			// Coordinates of node 1
			L1 = (double)(div_num - i) / div_num;
			L2 = (double)(i - j) / div_num;
			L3 = (double)j / div_num;
			sub_x1 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y1 = L1 * y1 + L2 * y2 + L3 * y3;
			// Coordinates of node 2
			L1 = (double)(div_num - i - 1) / div_num;
			L2 = (double)(i + 1 - j) / div_num;
			L3 = (double)j / div_num;
			sub_x2 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y2 = L1 * y1 + L2 * y2 + L3 * y3;
			// Coordinates of node 3
			L1 = (double)(div_num - i - 1) / div_num;
			L2 = (double)(i - j) / div_num;
			L3 = (double)(j + 1) / div_num;
			sub_x3 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y3 = L1 * y1 + L2 * y2 + L3 * y3;
			// Divided each small elements according to Hammer integration points
			/*
			pcl_param->x = 1.0 / 3.0 * sub_x1 + 1.0 / 3.0 * sub_x2 + 1.0 / 3.0 * sub_x3;
			pcl_param->y = 1.0 / 3.0 * sub_y1 + 1.0 / 3.0 * sub_y2 + 1.0 / 3.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area;
			addParticle(pcl_param, pcl_cm);
			*/
			// Integration point 1
			pcl_param->x = 2.0 / 3.0 * sub_x1 + 1.0 / 6.0 * sub_x2 + 1.0 / 6.0 * sub_x3;
			pcl_param->y = 2.0 / 3.0 * sub_y1 + 1.0 / 6.0 * sub_y2 + 1.0 / 6.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
			// Integration point 2
			pcl_param->x = 1.0 / 6.0 * sub_x1 + 2.0 / 3.0 * sub_x2 + 1.0 / 6.0 * sub_x3;
			pcl_param->y = 1.0 / 6.0 * sub_y1 + 2.0 / 3.0 * sub_y2 + 1.0 / 6.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
			// Integration point 3
			pcl_param->x = 1.0 / 6.0 * sub_x1 + 1.0 / 6.0 * sub_x2 + 2.0 / 3.0 * sub_x3;
			pcl_param->y = 1.0 / 6.0 * sub_y1 + 1.0 / 6.0 * sub_y2 + 2.0 / 3.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
		}

		for (j = 0; j < i; j++)
		{
			// Coordinates of node 1
			L1 = (double)(div_num - i - 1) / div_num;
			L2 = (double)(i - j) / div_num;
			L3 = (double)(j + 1) / div_num;
			sub_x1 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y1 = L1 * y1 + L2 * y2 + L3 * y3;
			// Coordinates of node 2
			L1 = (double)(div_num - i) / div_num;
			L2 = (double)(i - j - 1) / div_num;
			L3 = (double)(j + 1) / div_num;
			sub_x2 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y2 = L1 * y1 + L2 * y2 + L3 * y3;
			// Coordinates of node 3
			L1 = (double)(div_num - i) / div_num;
			L2 = (double)(i - j) / div_num;
			L3 = (double)j / div_num;
			sub_x3 = L1 * x1 + L2 * x2 + L3 * x3;
			sub_y3 = L1 * y1 + L2 * y2 + L3 * y3;
			// Divided each small elements according to Hammer integration points
			/*
			pcl_param->x = 1.0 / 3.0 * sub_x1 + 1.0 / 3.0 * sub_x2 + 1.0 / 3.0 * sub_x3;
			pcl_param->y = 1.0 / 3.0 * sub_y1 + 1.0 / 3.0 * sub_y2 + 1.0 / 3.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area;
			addParticle(pcl_param, pcl_cm);
			*/
			// Integration point 1
			pcl_param->x = 2.0 / 3.0 * sub_x1 + 1.0 / 6.0 * sub_x2 + 1.0 / 6.0 * sub_x3;
			pcl_param->y = 2.0 / 3.0 * sub_y1 + 1.0 / 6.0 * sub_y2 + 1.0 / 6.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
			// Integration point 2
			pcl_param->x = 1.0 / 6.0 * sub_x1 + 2.0 / 3.0 * sub_x2 + 1.0 / 6.0 * sub_x3;
			pcl_param->y = 1.0 / 6.0 * sub_y1 + 2.0 / 3.0 * sub_y2 + 1.0 / 6.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
			// Integration point 3
			pcl_param->x = 1.0 / 6.0 * sub_x1 + 1.0 / 6.0 * sub_x2 + 2.0 / 3.0 * sub_x3;
			pcl_param->y = 1.0 / 6.0 * sub_y1 + 1.0 / 6.0 * sub_y2 + 2.0 / 3.0 * sub_y3;
			pcl_param->mass = pcl_param->density * sub_elem_area / 3.0;
			addParticle(pcl_param, pcl_cm);
		}
	}
}