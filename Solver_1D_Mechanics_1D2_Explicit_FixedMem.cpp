#ifdef _DEBUG
#include <iostream>
#endif

#include "Solver_1D_Mechanics_1D2_Explicit_FixedMem.h"

#include "GaussIntegrationMacro.h"


Solver_1D_Mechanics_1D2_Explicit_FixedMem::Solver_1D_Mechanics_1D2_Explicit_FixedMem(
	const double time_step,
	Mesh_1D2 &mh,
	std::vector<ObjectByParticle_1D_Mechanics> &pcl_objs,
	OutputRequest &out) :
	mesh(mh), pcl_objects(pcl_objs),
	Solver(time_step, out),
	nodeVarMem(nullptr), particleVarMem(nullptr)
{
	nodeCoords = mesh.nodeCoords;
	nodes = mesh.nodes;
	nodeNum = mesh.nodeNum;
	elements = mesh.elements;
	elementNum = mesh.elementNum;

	object_num = pcl_objects.size();
}

Solver_1D_Mechanics_1D2_Explicit_FixedMem::
	~Solver_1D_Mechanics_1D2_Explicit_FixedMem()
{
	if (nodeVarMem) delete[] nodeVarMem;
	if (particleVarMem) delete[] particleVarMem;
}

void Solver_1D_Mechanics_1D2_Explicit_FixedMem::CFLTimeStep(void)
{
	//dt = 0.0; // need further work here...
}

// allocate space for all particle variables
int Solver_1D_Mechanics_1D2_Explicit_FixedMem::init(void)
{
	size_t i, j, k;
	size_t total_pcl_num, pcl_num;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	// Initialize calculation variables at nodes
	nodeVarMem = new NodeVar_1D_Mechanics[nodeNum * object_num];
	if (!nodeVarMem) throw std::exception("Fail to allocate nodal variables");

	for (i = 0; i < nodeNum; i++)
	{
		pNodeTmp = mesh.getNodeById(i + 1);
		pNodeVarTmp = nodeVarMem + i * object_num;
		// node and calculation variables point to each other
		pNodeTmp->nodeVar = pNodeVarTmp;
		for (j = 0; j < object_num; j++)
			pNodeVarTmp[j].node = pNodeTmp;
	}

	// Initialize calculation variables at particles
	total_pcl_num = 0;
	for (i = 0; i < object_num; i++)
		total_pcl_num += pcl_objects[i].particleNum;
	particleVarMem = new ParticleVar_1D_Mechanics[total_pcl_num];
	if (!particleVarMem) throw std::exception("Fail to allocate particle variables");

	k = 0;
	for (i = 0; i < object_num; i++)
	{
		curObject = &pcl_objects[i];
		pcl_num = curObject->particleNum;
		for (j = 0; j < pcl_num; j++)
		{
			pPclTmp = curObject->getParticleById(j + 1);
			pPclVarTmp = particleVarMem + k;
			++k;
			// particle and calculation variables point to each other
			pPclTmp->particleVar = pPclVarTmp;
			pPclVarTmp->particle = pPclTmp;
			// Initialize
			pPclTmp->isInMesh = true;
			pPclVarTmp->element = nullptr;
		}
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit_FixedMem::iteration(void)
{
	size_t i;

	for (i = 0; i < object_num; i++)
	{
		curObject = &pcl_objects[i];
		initStep_SingleObject();
		mapPhysicalProperyToNode_SingleObject();
		calInternalForce_SingleObject();
		calExternalForce_SingleObject();
		updatePhysicalPropertyAtNode_SingleObject();
		mapPhysicalPropertyToParticle_SingleObject();
		updatePhysicalPropertyAtParticle_SingleObject();
	}

	return 0;
}

// 1. Find in which elements where each particles lie;
// 3. Allocate space for node variables;
// 2. Calculate shape function and its one order derivative;
int Solver_1D_Mechanics_1D2_Explicit_FixedMem::initStep_SingleObject()
{
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	Element_1D2 *pElemTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;
	size_t i;
	size_t obj_index;
	size_t pcl_num_tmp;
	double xi;
	double x1, x2;
	double dN1_dxi, dN2_dxi;
	double dx_dxi, dxi_dx;

	obj_index = curObject->index - 1;
	
	// init all calculation variables
	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(nodes[i].nodeVar) + obj_index;
		pNodeVarTmp->needCal = false;
	}
	
	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			//std::cout << pPclVarTmp->particle->index << std::endl;
			
			pElemTmp = mesh.findInWhichElement(pPclTmp->x,
				static_cast<Element_1D2 *>(pPclVarTmp->element));
			if (pElemTmp) // if particles still lie in the mesh
			{
				pPclVarTmp->element = pElemTmp;
				mesh.calNaturalCoords(pElemTmp, pPclTmp->x, &xi);
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

				// node 1
				pNodeTmp = pElemTmp->node1;
				pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar) + obj_index;
				pPclVarTmp->nodeVar1 = pNodeVarTmp;
				//std::cout << pPclVarTmp->nodeVar1 << std::endl;
				// Initialize node calculation variables
				if (!(pNodeVarTmp->needCal))
				{
					pNodeVarTmp->needCal = true;
					pNodeVarTmp->mass = 0.0;
					pNodeVarTmp->momentum1 = 0.0;
					pNodeVarTmp->internalForce1 = 0.0;
					pNodeVarTmp->externalForce1 = 0.0;
					pNodeVarTmp->contactForce1 = 0.0;
				}

				// node 2
				pNodeTmp = pElemTmp->node2;
				pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar) + obj_index;
				pPclVarTmp->nodeVar2 = pNodeVarTmp;
				//std::cout << pPclVarTmp->nodeVar2 << std::endl;
				if (!(pNodeVarTmp->needCal))
				{
					pNodeVarTmp->needCal = true;
					pNodeVarTmp->mass = 0.0;
					pNodeVarTmp->momentum1 = 0.0;
					pNodeVarTmp->internalForce1 = 0.0;
					pNodeVarTmp->externalForce1 = 0.0;
					pNodeVarTmp->contactForce1 = 0.0;
				}
			}
			else
			{
				pPclTmp->isInMesh = false;
				pPclTmp->particleVar = nullptr;
			}
		}
	}

	return 0;
}


/*
This function is the first step of each calculation iteration, it:
1. Map mass and linear momentum from particles to nodes.
*/
int Solver_1D_Mechanics_1D2_Explicit_FixedMem::
	mapPhysicalProperyToNode_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		//std::cout << pPclTmp->index << std::endl;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			//std::cout << pNodeVarTmp << std::endl;
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N1;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N1;
			// node 2
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			//std::cout << pNodeVarTmp << std::endl;
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N2;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N2;
		}
	}

	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit_FixedMem::
	calInternalForce_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				pPclVarTmp->dN1_dx * pPclTmp->stress11;
			// node 2
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				pPclVarTmp->dN2_dx * pPclTmp->stress11;
		}
	}

	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit_FixedMem::
	calExternalForce_SingleObject(void)
{
	NodeVar_1D_Mechanics *pNodeVarTmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	MassForceBC *mfp_iter;
	double mass_force_tmp;
	SurfaceForceBC *sfp_iter;
	double surface_force_tmp;

	// ------------------- Boundary conditions specified on particles ------------------- 
	// body force per unit mass specified on particles 
	for (mfp_iter = curObject->massForceBCs_mem.get_first(); mfp_iter;
		 mfp_iter = curObject->massForceBCs_mem.get_next(mfp_iter))
	{
		pPclTmp = curObject->getParticleById(mfp_iter->index);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			if (mfp_iter->dof == DegreeOfFreedom::x)
			{
				mass_force_tmp = mfp_iter->massForce(t_cal);
				// node 1
				pNodeVarTmp = pPclVarTmp->nodeVar1;
				pNodeVarTmp->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N1;
				// node 2
				pNodeVarTmp = pPclVarTmp->nodeVar2;
				pNodeVarTmp->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N2;
			}
		}
	}

	//std::cout << "time: " << t_cal << std::endl;
	// surface force specified on particles at the surface of the objects 
	for (sfp_iter = curObject->surfaceForceBCs_mem.get_first(); sfp_iter;
		 sfp_iter = curObject->surfaceForceBCs_mem.get_next(sfp_iter))
	{
		pPclTmp = curObject->getParticleById(sfp_iter->index);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			if (sfp_iter->dof == DegreeOfFreedom::x)
			{
				surface_force_tmp = sfp_iter->surfaceForce(t_cal);
				//std::cout << "surfaceForce1: " << sfp_iter->surfaceForce1(t_cal) << std::endl;
				// node 1
				pNodeVarTmp = pPclVarTmp->nodeVar1;
				pNodeVarTmp->externalForce1 += surface_force_tmp * pPclVarTmp->N1;
				// node 2
				pNodeVarTmp = pPclVarTmp->nodeVar2;
				pNodeVarTmp->externalForce1 += surface_force_tmp * pPclVarTmp->N2;
			}
		}
	}

	return 0;
}


int Solver_1D_Mechanics_1D2_Explicit_FixedMem::
	updatePhysicalPropertyAtNode_SingleObject(void)
{
	size_t i;
	size_t obj_index;
	Node_1D *pNodeTmp;
	NodeVar_1D_Mechanics *pNodeVarTmp;
	AccelerationBC *a_iter;
	VelocityBC *v_iter;

	obj_index = curObject->index - 1;

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(nodes[i].nodeVar)	+ obj_index;
		if (pNodeVarTmp->needCal)
		{
			// calculate nodal force
			pNodeVarTmp->nodalForce1 = pNodeVarTmp->externalForce1
				- pNodeVarTmp->internalForce1 + pNodeVarTmp->contactForce1;
			// calculate increment of linear momentum
			pNodeVarTmp->dMomentum1 = pNodeVarTmp->nodalForce1 * t_increment_a;
		}
	}

	// Apply acceleration boundary conditions
	for (a_iter = mesh.accelerationBCs_mem.get_first(); a_iter;
		 a_iter = mesh.accelerationBCs_mem.get_next(a_iter))
	{
		pNodeTmp = mesh.getNodeById(a_iter->index);
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar) + obj_index;
		if (pNodeVarTmp->needCal)
		{
			if (a_iter->dof == DegreeOfFreedom::x)
			{
				pNodeVarTmp->dMomentum1 = a_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(nodes[i].nodeVar)	+ obj_index;
		if (pNodeVarTmp->needCal)
		{
			// integrate linear momentum
			pNodeVarTmp->momentum1 += pNodeVarTmp->dMomentum1;
		}
	}

	//std::cout << "time:" << t_cal << std::endl;
	// Apply velocity boundary conditions
	for (v_iter = mesh.velocityBCs_mem.get_first(); v_iter;
		 v_iter = mesh.velocityBCs_mem.get_next(v_iter))
	{
		pNodeTmp = mesh.getNodeById(v_iter->index);
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(pNodeTmp->nodeVar) + obj_index;
		if (pNodeVarTmp->needCal)
		{
			if (v_iter->dof == DegreeOfFreedom::x)
			{
				pNodeVarTmp->dMomentum1 = v_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
				pNodeVarTmp->momentum1 = v_iter->v(t_cal) * pNodeVarTmp->mass;
				//std::cout << "v1: " << v_iter->v1(t_cal) << std::endl;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Mechanics *>(nodes[i].nodeVar)	+ obj_index;
		if (pNodeVarTmp->needCal)
		{
			// calculate velocity
			pNodeVarTmp->v1 = pNodeVarTmp->momentum1 / pNodeVarTmp->mass;
			// calculate increment of displacement
			pNodeVarTmp->u = pNodeVarTmp->v1 * t_increment;
		}
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit_FixedMem::calContactForce(void)
{
	return 0;
}

// map momentum, displacement and strain
int Solver_1D_Mechanics_1D2_Explicit_FixedMem::mapPhysicalPropertyToParticle_SingleObject()
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	NodeVar_1D_Mechanics *pNode1VarTmp, *pNode2VarTmp;
	double N1Tmp, N2Tmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = static_cast<Particle_1D_Mechanics *>(curObject->particles) + i;
		
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);

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
	}

	return 0;
}

int Solver_1D_Mechanics_1D2_Explicit_FixedMem::updatePhysicalPropertyAtParticle_SingleObject()
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_1D_Mechanics *pPclTmp;
	ParticleVar_1D_Mechanics *pPclVarTmp;
	double dstrain[6], stress[6];
	double dstress[6], destrain[6], dpstrain[6];
	double F_det;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		//std::cout << pPclTmp->index << std::endl;

		if (pPclTmp->isInMesh) // only update element when in mesh
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Mechanics *>(pPclTmp->particleVar);
			
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
			//std::cout << (unsigned int)(pPclTmp->cm->getType()) << std::endl;
			pPclTmp->cm->integration_explicit(dstrain, stress, dstress, destrain, dpstrain);
			//std::cout << dstrain[0] << " " << dstress[0] << std::endl;

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

