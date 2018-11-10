#include "stdafx.h"
#include "Solver_1D_Mechanics_1D2_Explicit.h"

#include "GaussIntegrationMacro.h"

Solver_1D_Mechanics_1D2_Explicit::Solver_1D_Mechanics_1D2_Explicit(
	Data_1D_Mechanics_1D2 &da, OutputRequest &out) :
	data(da), mesh(da.mesh), objects(da.objectByParticle),
	Solver(da, out, "Solver_1D_Mechanics_1D2_Explicit")
{
	size_t i;
	nodeCoords = mesh.nodeCoords;
	nodes = mesh.nodes;
	nodeNum = mesh.nodeNum;
	elements = mesh.elements;
	elementNum = mesh.elementNum;
	size_t pcl_num;

	nodeVarPool.setSlotSize(nodeNum / 4); // according to "experience"
	node_contact.reserve(nodeNum / 10); // according to "experience"

	object_num = objects.size();
	pcl_num = 0;
	for (i = 0; i < object_num; i++)
	{
		pcl_num += objects[i].getParticleNum();
	}
	particleVarPool.setSlotSize(pcl_num);
	
	/*
	reset_stack();
	for (i = 0; i < nodeNum; i++)
		nodes[i].reset_stack();

	obj_num_tmp = objects.size();
	for (i = 0; i < obj_num_tmp; i++)
	{
		objects[i].reset_list();
	}*/
}
	
Solver_1D_Mechanics_1D2_Explicit::~Solver_1D_Mechanics_1D2_Explicit()
{
	particleVarPool.clearMem();
	nodeVarPool.clearMem();
}

void Solver_1D_Mechanics_1D2_Explicit::CFLTimeStep(void)
{
	//dt = 0.0; // need further work here...
}

// allocate space for all particle variables
int Solver_1D_Mechanics_1D2_Explicit::init(void)
{
	size_t i, j;
	size_t pcl_num_tmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;

	for (i = 0; i < object_num; i++)
	{
		curObject = &objects[i];
		pcl_num_tmp = curObject->particleNum;
		for (j = 0; j < pcl_num_tmp; j++)
		{
			pPclTmp = &(curObject->particles[j]);
			
			pPclVarTmp = particleVarPool.alloc();
			pPclTmp->particleVar = pPclVarTmp;
			pPclVarTmp->particle = pPclTmp;
			
			curObject->add_list(pPclVarTmp);
			
			pPclTmp->isInMesh = true;
			pPclVarTmp->element = nullptr;
		}
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit::iteration(void)
{
	initStep_SingleObject();
	mapPhysicalProperyToNode_SingleObject();
	calInternalForce_SingleObject();
	calExternalForce_SingleObject();
	updatePhysicalPropertyAtNode_noContact();
	mapPhysicalPropertyToParticle_SingleObject();
	updatePhysicalPropertyAtParticle_SingleObject();

	return 0;
}

// 1. Find in which elements where each particles lie;
// 3. Allocate space for node variables;
// 2. Calculate shape function and its one order derivative;
int Solver_1D_Mechanics_1D2_Explicit::initStep_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	Element_1D2 *pElemTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	double xi;
	double x1, x2;
	double dN1_dxi, dN2_dxi;
	double dx_dxi, dxi_dx;

	curObject->start_list();
	while(pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list()))
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		pElemTmp = mesh.findInWhichElement(pPclTmp->x,
			static_cast<Element_1D2 *>(pPclVarTmp->element));
		if (pElemTmp) // if particles still lie in the mesh
		{
			pPclVarTmp->element = pElemTmp;
			mesh.calNaturalCoords(pElemTmp, pPclTmp->x, &xi);
			//std::cout << "Natural coord: " << xi << " ";
			// Calculate shape functions
			pPclVarTmp->N1 = Mesh_1D2::N1(xi);
			pPclVarTmp->N2 = Mesh_1D2::N2(xi);
			// Calculate one order derivative of shape functions
			x1 = nodeCoords[pElemTmp->xIndex];
			x2 = nodeCoords[pElemTmp->xIndex + 1];
			dN1_dxi = Mesh_1D2::dN1_dxi();
			dN2_dxi = Mesh_1D2::dN2_dxi();
			dx_dxi = x1 * dN1_dxi + x2 * dN2_dxi;
			dxi_dx = 1.0 / dx_dxi;
			pPclVarTmp->dN1_dx = dN1_dxi * dxi_dx;
			pPclVarTmp->dN2_dx = dN2_dxi * dxi_dx;

			// Allocate space for node variables
			pNodeTmp = pElemTmp->node1;
			pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->get_top_stack());
			if (!pNodeVarTmp || pNodeVarTmp->objectIndex != curObject->index)
			{
				if (pNodeTmp->get_num_stack() == 1)
					node_contact.push_back(pNodeTmp->index);

				pNodeVarTmp = nodeVarPool.alloc();
				pNodeVarTmp->node = pNodeTmp;
				pNodeTmp->push_stack(pNodeVarTmp);
				push_stack(pNodeVarTmp);

				pNodeVarTmp->objectIndex = curObject->index;
				pNodeVarTmp->mass = 0.0;
				pNodeVarTmp->momentum1 = 0.0;
				pNodeVarTmp->internalForce1 = 0.0;
				pNodeVarTmp->externalForce1 = 0.0;

				// This node may be possible to have collision
				if (pNodeTmp->get_num_stack() == 2)
					node_contact.push_back(pNodeTmp->index);
				pNodeVarTmp->contactForce1 = 0.0;
			}
			pPclVarTmp->nodeVar1 = pNodeVarTmp;

			pNodeTmp = pElemTmp->node2;
			pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->get_top_stack());
			if (!pNodeVarTmp || pNodeVarTmp->objectIndex != curObject->index)
			{
				pNodeVarTmp = nodeVarPool.alloc();
				pNodeVarTmp->node = pNodeTmp;
				pNodeTmp->push_stack(pNodeVarTmp);
				push_stack(pNodeVarTmp);

				pNodeVarTmp->objectIndex = curObject->index;
				pNodeVarTmp->mass = 0.0;
				pNodeVarTmp->momentum1 = 0.0;
				pNodeVarTmp->internalForce1 = 0.0;
				pNodeVarTmp->externalForce1 = 0.0;
				
				// This node may be possible to have collision
				if (pNodeTmp->get_num_stack() == 2)
					node_contact.push_back(pNodeTmp->index);
				pNodeVarTmp->contactForce1 = 0.0;
				pNodeVarTmp->normal1 = 0.0;
			}
			pPclVarTmp->nodeVar2 = pNodeVarTmp;

			curObject->next_list();
		}
		else
		{
			pPclTmp->isInMesh = false;
			curObject->del_cur_list();
		}
	}

	return 0;
}


/*
This function is the first step of each calculation iteration, it:
1. Map mass and linear momentum from particles to nodes.
*/
int Solver_1D_Mechanics_1D2_Explicit::mapPhysicalProperyToNode_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	// loop for all particles in the object
	for (curObject->start2_list();
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list());
		curObject->next2_list())
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		// node 1
		pNodeVarTmp = pPclVarTmp->nodeVar1;
		pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N1;
		pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N1;
		// node 2
		pNodeVarTmp = pPclVarTmp->nodeVar2;
		pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N2;
		pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N2;
	}

	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit::calInternalForce_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	// loop for all particles in the object
	for (curObject->start2_list();
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list());
		curObject->next2_list())
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		// node 1
		pNodeVarTmp = pPclVarTmp->nodeVar1;
		pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
			pPclVarTmp->dN1_dx * pPclTmp->stress11;
		// node 2
		pNodeVarTmp = pPclVarTmp->nodeVar2;
		pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
			pPclVarTmp->dN2_dx * pPclTmp->stress11;
	}
	
	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit::calExternalForce_SingleObject()
{
	size_t i;
	NodeVar_1D_Mechanics *pNodeVarTmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	size_t massForceOnParticleNum;
	MassForceOnParticle_1D_Mechanics *massForceOnParticle;
	size_t surfaceForceOnParticleNum;
	SurfaceForceOnParticle_1D_Mechanics *surfaceForceOnParticle;

	// ------------------- Boundary conditions specified on particles ------------------- 
	// body force per unit mass specified on particles 
	massForceOnParticleNum = curObject->massForceOnParticleNum;
	massForceOnParticle = curObject->massForceOnParticle;
	for (i = 0; i < massForceOnParticleNum; ++i)
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(massForceOnParticle[i].particle);
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
		// node 1
		pNodeVarTmp = pPclVarTmp->nodeVar1;
		pNodeVarTmp->externalForce1 += pPclTmp->mass * massForceOnParticle[i].bodyForce1
			* pPclVarTmp->N1;
		// node 2
		pNodeVarTmp = pPclVarTmp->nodeVar2;
		pNodeVarTmp->externalForce1 += pPclTmp->mass * massForceOnParticle[i].bodyForce1
			* pPclVarTmp->N2;
	}

	// surface force specified on particles at the surface of the objects 
	surfaceForceOnParticleNum = curObject->surfaceForceOnParticleNum;
	surfaceForceOnParticle = curObject->surfaceForceOnParticle;
	for (i = 0; i < surfaceForceOnParticleNum; i++)
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(surfaceForceOnParticle[i].particle);
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
		// node 1
		pNodeVarTmp = pPclVarTmp->nodeVar1;
		pNodeVarTmp->externalForce1 += surfaceForceOnParticle[i].surfaceForce1 * pPclVarTmp->N1;
		// node 2
		pNodeVarTmp = pPclVarTmp->nodeVar2;
		pNodeVarTmp->externalForce1 += surfaceForceOnParticle[i].surfaceForce1 * pPclVarTmp->N2;
	}

	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit::updatePhysicalPropertyAtNode_noContact()
{
	size_t i;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;
	size_t accelerationBCNum;
	AccelerationBC_1D *accelerationBCs;
	size_t velocityBCNum;
	VelocityBC_1D *velocityBCs;

	for (start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)); 
		pNodeVarTmp; next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
	{
		// calculate nodal force
		pNodeVarTmp->nodalForce1 = pNodeVarTmp->externalForce1 - pNodeVarTmp->internalForce1;
		// calculate increment of linear momentum
		pNodeVarTmp->dMomentum1 = pNodeVarTmp->nodalForce1 * t_increment_a;
	}

	// Apply acceleration boundary conditions
	accelerationBCNum = data.accelerationBCNum;
	accelerationBCs = data.accelerationBCs;
	for (i = 0; i < accelerationBCNum; i++)
	{
		pNodeTmp = static_cast<Node_1D *>(accelerationBCs[i].node);
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar);
		pNodeVarTmp->dMomentum1 = accelerationBCs[i].a1 * pNodeVarTmp->mass * t_increment_a;
	}

	for (start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp));
		pNodeVarTmp; next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
	{
		// integrate linear momentum
		pNodeVarTmp->momentum1 += pNodeVarTmp->dMomentum1;
	}

	// Apply velocity boundary conditions
	velocityBCNum = data.velocityBCNum;
	velocityBCs = data.velocityBCs;
	for (i = 0; i < velocityBCNum; i++)
	{
		pNodeTmp = static_cast<Node_1D *>(velocityBCs[i].node);
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar);
		// assume acceleration = 0 at velocity boundary conditions
		pNodeVarTmp->dMomentum1 = 0.0;
		pNodeVarTmp->momentum1 = velocityBCs[i].v1 * pNodeVarTmp->mass;
	}

	for (start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp));
		pNodeVarTmp; next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
	{
		// calculate velocity
		pNodeVarTmp->v1 = pNodeVarTmp->momentum1 / pNodeVarTmp->mass;
		// calculate increment of displacement
		pNodeVarTmp->u = pNodeVarTmp->v1 * t_increment;
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit::calContactForce(void)
{
	size_t i, j;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;
	double mass_sum_tmp, momentum_sum_tmp;
	size_t node_contact_num;
	std::vector<bool> is_contact;
	size_t contact_num;

	// calculate normal vector
	// for each object
	for (curObject->start2_list();
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list());
		curObject->next2_list())
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		// node 1
		pNodeVarTmp = pPclVarTmp->nodeVar1;
		pNodeTmp = static_cast<Node_1D *>(pNodeVarTmp->node);
		if (pNodeTmp->get_num_stack() > 1)
			pNodeVarTmp->normal1 += pPclTmp->mass * pPclVarTmp->dN1_dx;
		// node 2
		pNodeVarTmp = pPclVarTmp->nodeVar2;
		pNodeTmp = static_cast<Node_1D *>(pNodeVarTmp->node);
		if (pNodeTmp->get_num_stack() > 1)
			pNodeVarTmp->normal1 += pPclTmp->mass * pPclVarTmp->dN2_dx;
	}
	node_contact_num = node_contact.size();
	for (i = 0; i < node_contact_num; i++)
	{
		pNodeTmp = mesh.getNodeById(node_contact[i]);
		for (pNodeTmp->start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)); pNodeVarTmp;
			pNodeTmp->next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
		{
			// normalized normal vector;
			pNodeVarTmp->normal1 = pNodeVarTmp->normal1 >= 0 ? 1.0 : -1.0;
		}
	}

	node_contact_num = node_contact.size();
	is_contact.reserve(5);
	for (i = 0; i < node_contact_num; i++)
	{
		pNodeTmp = mesh.getNodeById(node_contact[i]);
		mass_sum_tmp = 0.0;
		momentum_sum_tmp = 0.0;
		// initial node var
		for (pNodeTmp->start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)); pNodeVarTmp;
			pNodeTmp->next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
		{
			mass_sum_tmp += pNodeVarTmp->mass;
			momentum_sum_tmp += pNodeVarTmp->momentum1;
		}
		// calculate center of mass velocity
		pNodeTmp->v1_cm = momentum_sum_tmp / mass_sum_tmp;

		is_contact.clear();
		contact_num = pNodeTmp->get_num_stack();
		for (pNodeTmp->start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)); pNodeVarTmp;
			pNodeTmp->next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
		{
			// detect collision
			if ((pNodeVarTmp->v1 - pNodeTmp->v1_cm) * pNodeVarTmp->normal1 > 0)
			{
				is_contact.push_back(true);
			}
			else
			{
				--contact_num;
				is_contact.push_back(false);
				mass_sum_tmp -= pNodeVarTmp->mass;
				momentum_sum_tmp -= pNodeVarTmp->momentum1;
			}
		}

		if (contact_num)
		{
			pNodeTmp->v1_cm = momentum_sum_tmp / mass_sum_tmp;
			j = 0;
			for (pNodeTmp->start_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)); pNodeVarTmp;
				pNodeTmp->next_stack(reinterpret_cast<NodeVar **>(&pNodeVarTmp)))
			{
				if (is_contact[j])
					pNodeVarTmp->contactForce1 = (pNodeTmp->v1_cm - pNodeVarTmp->v1)
						* pNodeVarTmp->mass / t_increment_a;
				++j;
			}
		}

	}

	return 0;
}

// map momentum, displacement and strain
int Solver_1D_Mechanics_1D2_Explicit::mapPhysicalPropertyToParticle_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNode1VarTmp, *pNode2VarTmp;
	double N1Tmp, N2Tmp;

	for (curObject->start2_list();
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list());
		curObject->next2_list())
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		// --------------------------- node 1 ------------------------------
		pNode1VarTmp = static_cast<NodeVar_1D_Mechanics *>(pPclVarTmp->nodeVar1);
		N1Tmp = pPclVarTmp->N1;
		// --------------------------- node 2 ------------------------------
		pNode2VarTmp = static_cast<NodeVar_1D_Mechanics *>(pPclVarTmp->nodeVar2);
		N2Tmp = pPclVarTmp->N2;

		// update linear momentum of particle
		pPclTmp->momentum1 += pPclTmp->mass *
				(pNode1VarTmp->dMomentum1 / pNode1VarTmp->mass * N1Tmp
			+ pNode2VarTmp->dMomentum1 / pNode2VarTmp->mass * N2Tmp);
		// update location of particle
		pPclTmp->x += pNode1VarTmp->u * N1Tmp + pNode2VarTmp->u * N2Tmp;
		// update strain increment of particle
		pPclVarTmp->dstrain11 =
			pNode1VarTmp->u * pPclVarTmp->dN1_dx + pNode2VarTmp->u * pPclVarTmp->dN2_dx;
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit::updatePhysicalPropertyAtParticle_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	double dstrain[6], stress[6];
	double dstress[6], destrain[6], dpstrain[6];
	double F_det;

	for (curObject->start2_list();
		pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(curObject->get_cur_list());
		curObject->next2_list())
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(pPclVarTmp->particle);
		if (pPclTmp->isInMesh) // only update element when in mesh
		{
			// update density
			F_det = 1.0 + pPclVarTmp->dstrain11;
			pPclTmp->density /= F_det;

			// integrate constitute
			dstrain[0] = pPclVarTmp->dstrain11;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			stress[0] = pPclTmp->stress11;
			stress[1] = 0.0;
			stress[2] = 0.0;
			stress[3] = 0.0;
			stress[4] = 0.0;
			stress[5] = 0.0;
			
			pPclTmp->cm->integration_explicit(dstrain, stress, dstress, destrain, dpstrain);
			
			//update stress
			pPclTmp->stress11 += dstress[0];
			// update strain (also assume that the strain increment is Jaumann rate)
			pPclTmp->strain11 += pPclVarTmp->dstrain11;
			pPclTmp->estrain11 += destrain[0];
			pPclTmp->pstrain11 += dpstrain[0];
		}
	}

	return 0;
}
