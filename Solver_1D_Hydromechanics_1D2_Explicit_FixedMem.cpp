#ifdef _DEBUG
#include <iostream>
#endif

#include "Solver_1D_Hydromechanics_1D2_Explicit_FixedMem.h"

#include "GaussIntegrationMacro.h"


Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(
	double time_step,
	Mesh_1D2 &mh,
	std::vector<ObjectByParticle_1D_Hydromechanics> &pcl_objs,
	OutputRequest &out,
	const char *na) :
	Solver(time_step, out, na),
	mesh(mh), pcl_objects(pcl_objs),
	nodeVarMem(nullptr), particleVarMem(nullptr)
{
	nodeCoords = mesh.nodeCoords;
	nodes = mesh.nodes;
	nodeNum = mesh.nodeNum;
	elements = mesh.elements;
	elementNum = mesh.elementNum;

	object_num = pcl_objects.size();
}


Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(
		double time_step,
		Solver_1D_Hydromechanics_1D2_Explicit_FixedMem &prev_solver,
		const char *na) :
	Solver(time_step, prev_solver, na),
	mesh(prev_solver.mesh), pcl_objects(prev_solver.pcl_objects),
	nodeVarMem(prev_solver.nodeVarMem),
	particleVarMem(prev_solver.particleVarMem)
{
	nodeCoords = mesh.nodeCoords;
	nodes = mesh.nodes;
	nodeNum = mesh.nodeNum;
	elements = mesh.elements;
	elementNum = mesh.elementNum;

	object_num = pcl_objects.size();

	prev_solver.nodeVarMem = nullptr;
	prev_solver.particleVarMem = nullptr;
}


Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	~Solver_1D_Hydromechanics_1D2_Explicit_FixedMem()
{
	if (nodeVarMem) delete[] nodeVarMem;
	if (particleVarMem) delete[] particleVarMem;
}


void Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::CFLTimeStep(void)
{
	//dt = 0.0; // need further work here...
}


int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::init(void)
{
	size_t i, j, k;
	size_t total_pcl_num, pcl_num;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;

	// if is the first step
	if (index == 1)
	{
		// Initialize calculation variables at nodes
		nodeVarMem = new NodeVar_1D_Hydromechanics[nodeNum * object_num];
		if (!nodeVarMem) throw std::exception("Fail to allocate nodal variables.");

		for (i = 0; i < nodeNum; i++)
		{
			pNodeTmp = mesh.nodes + i;
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
		particleVarMem = new ParticleVar_1D_Hydromechanics[total_pcl_num];
		if (!particleVarMem) throw std::exception("Fail to allocate particle variables");

		k = 0;
		for (i = 0; i < object_num; ++i)
		{
			curObject = &pcl_objects[i];
			pcl_num = curObject->particleNum;
			for (j = 0; j < pcl_num; ++j)
			{
				pPclTmp = curObject->particles + j;
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
	}

	return 0;
}


int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::iteration(void)
{
	size_t i;

	for (i = 0; i < object_num; i++)
	{
		curObject = &pcl_objects[i];
		curObjectId = i;
		initStep_SingleObject();
		mapPhysicalProperyToNode_SingleObject();
		calInternalForce_SingleObject();
		calExternalForce_SingleObject();
		updatePhysicalPropertyAtNode_SingleObject();
		calFluidPhase();
		calContactForce();
		mapPhysicalPropertyToParticle_SingleObject();
		updatePhysicalPropertyAtParticle_SingleObject();
	}
	
	return 0;
}

// 1. Find in which elements where each particles lie;
// 3. Allocate space for node variables;
// 2. Calculate shape function and its one order derivative;
int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	initStep_SingleObject(void)
{
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	Element_1D2 *pElemTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;
	size_t i;
	size_t pcl_num_tmp;
	double xi;
	double x1, x2;
	double dN1_dxi, dN2_dxi;
	double dx_dxi, dxi_dx;

	// init all calculation variables
	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar) + curObjectId;
		pNodeVarTmp->needCal = false;
	}

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			//std::cout << pPclVarTmp->particle->index << std::endl;

			pElemTmp = mesh.findInWhichElement(pPclTmp->x,
				static_cast<Element_1D2 *>(pPclVarTmp->element));
			if (pElemTmp) // if particles still lie in the mesh
			{
				// Initialize particle calculation variables
				pPclVarTmp->element = pElemTmp;
				mesh.calNaturalCoords(pElemTmp, pPclTmp->x, &xi);
				pPclVarTmp->N1 = Mesh_1D2::N1(xi);
				pPclVarTmp->N2 = Mesh_1D2::N2(xi);
				//std::cout << pPclVarTmp->N1 << " " << pPclVarTmp->N2 << std::endl;
				// Calculate one order derivative of shape functions
				x1 = nodeCoords[pElemTmp->xIndex];
				x2 = nodeCoords[pElemTmp->xIndex + 1];
				dN1_dxi = Mesh_1D2::dN1_dxi();
				dN2_dxi = Mesh_1D2::dN2_dxi();
				dx_dxi = x1 * dN1_dxi + x2 * dN2_dxi;
				dxi_dx = 1.0 / dx_dxi;
				pPclVarTmp->dN1_dx = dN1_dxi * dxi_dx;
				pPclVarTmp->dN2_dx = dN2_dxi * dxi_dx;

				pPclVarTmp->avgdensity_s = (1.0 - pPclTmp->n) * pPclTmp->density_s;
				pPclVarTmp->avgdensity_f = pPclTmp->n * pPclTmp->density_f;
				pPclVarTmp->volume = pPclTmp->mass_s / pPclVarTmp->avgdensity_s;

				/* ---------------------------------------------------------------
				Initialize nodal variables
				---------------------------------------------------------------- */
#define			INIT_NODEVAR(pNodeVar) \
				do { \
					if (!(pNodeVar->needCal)) \
					{ \
						pNodeVar->needCal = true; \
						pNodeVar->momentum1_m = 0.0; \
						pNodeVar->internalForce1_m = 0.0; \
						pNodeVar->externalForce1_m = 0.0; \
						pNodeVar->contactForce1_m = 0.0; \
						pNodeVar->volume = 0.0; \
						pNodeVar->mass_s = 0.0; \
						pNodeVar->mass_f = 0.0; \
						pNodeVar->internalForce1_f = 0.0; \
						pNodeVar->externalForce1_f = 0.0; \
					} \
				} while (0)

				// node 1
				pNodeTmp = pElemTmp->node1;
				pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar1 = pNodeVarTmp;

				// node 2
				pNodeTmp = pElemTmp->node2;
				pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar2 = pNodeVarTmp;
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
int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	mapPhysicalProperyToNode_SingleObject(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			pNodeVarTmp->momentum1_m += pPclVarTmp->N1 *pPclTmp->momentum1_m;
			pNodeVarTmp->mass_s += pPclVarTmp->N1 * pPclTmp->mass_s;
			pNodeVarTmp->mass_f += pPclVarTmp->N1 *pPclTmp->mass_f;
			pNodeVarTmp->volume += pPclVarTmp->N1 * pPclVarTmp->volume;
			// node 2
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			pNodeVarTmp->momentum1_m += pPclVarTmp->N2 * pPclTmp->momentum1_m;
			pNodeVarTmp->mass_s += pPclVarTmp->N2 * pPclTmp->mass_s;
			pNodeVarTmp->mass_f += pPclVarTmp->N2 * pPclTmp->mass_f;
			pNodeVarTmp->volume += pPclVarTmp->N2 * pPclVarTmp->volume;
		}
	}

	return 0;
}


int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	calInternalForce_SingleObject(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		if (pPclTmp->isInMesh)
		{
			//std::cout << "p_id: " << pPclTmp->index;
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			// ------------------- Node 1 ----------------------
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			pNodeVarTmp->internalForce1_m += pPclVarTmp->dN1_dx * (pPclTmp->estress11 - pPclTmp->p) * pPclVarTmp->volume;
			pNodeVarTmp->internalForce1_f += pPclVarTmp->dN1_dx * pPclTmp->k * -pPclTmp->p * pPclVarTmp->volume;
			// ------------------- Node 2 ----------------------
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			pNodeVarTmp->internalForce1_m += pPclVarTmp->dN2_dx * (pPclTmp->estress11 - pPclTmp->p) * pPclVarTmp->volume;
			pNodeVarTmp->internalForce1_f += pPclVarTmp->dN2_dx * pPclTmp->k * -pPclTmp->p * pPclVarTmp->volume;
			//std::cout << "dN2_dx: " << pPclVarTmp->dN2_dx << " p: " << pPclTmp->p << std::endl;
			//std::cout << " p: " << pPclTmp->p << std::endl;
		}
	}

	return 0;
}

int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	calExternalForce_SingleObject(void)
{
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	NodeVar_1D_Hydromechanics *pNode1VarTmp, *pNode2VarTmp;
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
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			pNode2VarTmp = pPclVarTmp->nodeVar2;
			mass_force_tmp = mfp_iter->massForce(t_cal);
			switch (mfp_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNode1VarTmp->externalForce1_m += pPclVarTmp->N1 * mass_force_tmp * pPclTmp->mass_s;
				pNode2VarTmp->externalForce1_m += pPclVarTmp->N2 * mass_force_tmp * pPclTmp->mass_s;
				break;
			case DegreeOfFreedom::x_f:
				pNode1VarTmp->externalForce1_m += pPclVarTmp->N1 * mass_force_tmp * pPclTmp->mass_f;
				pNode2VarTmp->externalForce1_m += pPclVarTmp->N2 * mass_force_tmp * pPclTmp->mass_f;
				pNode1VarTmp->externalForce1_f += pPclVarTmp->N1 * pPclTmp->k * mass_force_tmp * pPclTmp->density_f * pPclVarTmp->volume;
				pNode2VarTmp->externalForce1_f += pPclVarTmp->N2 * pPclTmp->k * mass_force_tmp * pPclTmp->density_f * pPclVarTmp->volume;
				break;
			default:
				break;
			}
		}
	}

	// surface force specified on particles at the surface of the objects 
	for (sfp_iter = curObject->surfaceForceBCs_mem.get_first(); sfp_iter;
		sfp_iter = curObject->surfaceForceBCs_mem.get_next(sfp_iter))
	{
		pPclTmp = curObject->getParticleById(sfp_iter->index);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			pNode2VarTmp = pPclVarTmp->nodeVar2;
			surface_force_tmp = sfp_iter->surfaceForce(t_cal);
			switch (sfp_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNode1VarTmp->externalForce1_m += pPclVarTmp->N1 * surface_force_tmp;
				pNode2VarTmp->externalForce1_m += pPclVarTmp->N2 * surface_force_tmp;
				break;
			case DegreeOfFreedom::x_f:
				pNode1VarTmp->externalForce1_m += pPclVarTmp->N1 * surface_force_tmp;
				pNode2VarTmp->externalForce1_m += pPclVarTmp->N2 * surface_force_tmp;
				pNode1VarTmp->externalForce1_f += pPclVarTmp->N1 * pPclTmp->k * surface_force_tmp;
				pNode2VarTmp->externalForce1_f += pPclVarTmp->N2 * pPclTmp->k * surface_force_tmp;
				break;
			default:
				break;
			}
		}
	}

	return 0;
}

int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	updatePhysicalPropertyAtNode_SingleObject(void)
{
	size_t i;
	Node_1D *pNodeTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;
	AccelerationBC *a_iter;
	VelocityBC *v_iter;

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			// calculate nodal force
			pNodeVarTmp->nodalForce1_m = pNodeVarTmp->externalForce1_m
				- pNodeVarTmp->internalForce1_m + pNodeVarTmp->contactForce1_m;
			
			//std::cout << "n_id: " << nodes[i].index
			//	<< " ef: " << pNodeVarTmp->externalForce1_m
			//	<< " if: " << pNodeVarTmp->internalForce1_m
			//	<< " cf: " << pNodeVarTmp->contactForce1_m
			//	<< " nf: " << pNodeVarTmp->nodalForce1_m << std::endl;

			// calculate increment of linear momentum
			pNodeVarTmp->dMomentum1_m = pNodeVarTmp->nodalForce1_m * t_increment_a;
			pNodeVarTmp->a1_s = pNodeVarTmp->nodalForce1_m / (pNodeVarTmp->mass_s + pNodeVarTmp->mass_f);
		
			//std::cout << "n_id: " << nodes[i].index
			//	<< " dm1: " << pNodeVarTmp->dMomentum1_m << std::endl;
		}
	}
	
	// Apply acceleration boundary conditions
	for (a_iter = mesh.accelerationBCs_mem.get_first(); a_iter;
		a_iter = mesh.accelerationBCs_mem.get_next(a_iter))
	{
		pNodeTmp = mesh.getNodeById(a_iter->index);
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(pNodeTmp->nodeVar);
		if (pNodeVarTmp->needCal)
		{
			switch (a_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNodeVarTmp->a1_s = a_iter->a(t_cal);
				pNodeVarTmp->dMomentum1_m = (pNodeVarTmp->mass_s + pNodeVarTmp->mass_f) * pNodeVarTmp->a1_s * t_increment_a;
				break;
			default:
				break;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar);
		if (pNodeVarTmp->needCal)
		{
			// update linear momentum
			pNodeVarTmp->momentum1_m += pNodeVarTmp->dMomentum1_m;
			//std::cout << "n_id: " << nodes[i].index
			//	<< " dm1: " << pNodeVarTmp->dMomentum1_m
			//	<< " m1: " << pNodeVarTmp->momentum1_m << std::endl;
		}
	}

	// Apply velocity boundary conditions
	for (v_iter = mesh.velocityBCs_mem.get_first(); v_iter;
		 v_iter = mesh.velocityBCs_mem.get_next(v_iter))
	{
		pNodeTmp = mesh.getNodeById(v_iter->index);
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(pNodeTmp->nodeVar);
		if (pNodeVarTmp->needCal)
		{
			switch (v_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNodeVarTmp->a1_s = v_iter->a(t_cal);
				pNodeVarTmp->dMomentum1_m = (pNodeVarTmp->mass_s + pNodeVarTmp->mass_f) * pNodeVarTmp->a1_s * t_increment_a;
				pNodeVarTmp->momentum1_m  = (pNodeVarTmp->mass_s + pNodeVarTmp->mass_f) * v_iter->v(t_cal);
				break;
			default:
				break;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar);
		if (pNodeVarTmp->needCal)
		{
			// calculate velocity
			pNodeVarTmp->v1_s = pNodeVarTmp->momentum1_m / (pNodeVarTmp->mass_s + pNodeVarTmp->mass_f);
			// calculate increment of displacement
			pNodeVarTmp->u1_s = pNodeVarTmp->v1_s * t_increment;
			//std::cout << "n_id: " << nodes[i].index 
			//	<< " v1_s: " << pNodeVarTmp->v1_s << std::endl;
		}
	}

	return 0;
}

int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	calFluidPhase(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	Node_1D *pNodeTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;
	NodeVar_1D_Hydromechanics *pNode1VarTmp, *pNode2VarTmp;
	VelocityBC *v_iter;

	// map acceleration back to particle
	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			pNode2VarTmp = pPclVarTmp->nodeVar2;

			// map acceleration back to particle
			pPclVarTmp->a1_s = pNode1VarTmp->a1_s * pPclVarTmp->N1 + pNode2VarTmp->a1_s * pPclVarTmp->N2;
			//std::cout << "a1_s: " << pNode1VarTmp->a1_s << " " << pNode2VarTmp->a1_s << std::endl;
			pNode1VarTmp->externalForce1_f += pPclVarTmp->N1 * pPclTmp->k * -pPclVarTmp->a1_s * pPclTmp->density_f * pPclVarTmp->volume;
			pNode2VarTmp->externalForce1_f += pPclVarTmp->N2 * pPclTmp->k * -pPclVarTmp->a1_s * pPclTmp->density_f * pPclVarTmp->volume;
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar);
		if (pNodeVarTmp->needCal)
		{
			pNodeVarTmp->w1 = (pNodeVarTmp->externalForce1_f - pNodeVarTmp->internalForce1_f) / pNodeVarTmp->volume;
			pNodeVarTmp->u1_f = pNodeVarTmp->w1 * t_increment;
			//std::cout << "n_id: " << nodes[i].index << " w1: " << pNodeVarTmp->w1
			//	<< " ef_f: " << pNodeVarTmp->externalForce1_f << " if_f: " << pNodeVarTmp->internalForce1_f << std::endl;
		}
	}

	// Apply velocity boundary conditions
	for (v_iter = mesh.velocityBCs_mem.get_first(); v_iter;
		v_iter = mesh.velocityBCs_mem.get_next(v_iter))
	{
		pNodeTmp = mesh.getNodeById(v_iter->index);
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(pNodeTmp->nodeVar);
		if (pNodeVarTmp->needCal)
		{
			switch (v_iter->dof)
			{
			case DegreeOfFreedom::x_f:
				pNodeVarTmp->w1 = v_iter->v(t_cal);
				pNodeVarTmp->u1_f = pNodeVarTmp->w1 * t_increment;
				break;
			default:
				break;
			}
		}
	}
	
	/*
	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar);
		if (pNodeVarTmp->needCal)
		{
			std::cout << "n_id: " << nodes[i].index
				<< " u_s: " << pNodeVarTmp->u1_s
				<< " u_f: " << pNodeVarTmp->u1_f << std::endl;
		}
	}*/

	return 0;
}

int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	calContactForce(void)
{
	return 0;
}

// map momentum, displacement and strain
int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	mapPhysicalPropertyToParticle_SingleObject(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	NodeVar_1D_Hydromechanics *pNode1VarTmp, *pNode2VarTmp;
	double N1Tmp, N2Tmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			N1Tmp = pPclVarTmp->N1;
			N2Tmp = pPclVarTmp->N2;
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			pNode2VarTmp = pPclVarTmp->nodeVar2;

			// map acceleration back to particle
			//pPclVarTmp->a1_s = pNode1VarTmp->a1_s * N1Tmp + pNode2VarTmp->a1_s * N2Tmp;
			// update linear momentum of soil - fluid mixture
			pPclTmp->momentum1_m += (pPclTmp->mass_s + pPclTmp->mass_f) *
				(pNode1VarTmp->dMomentum1_m / (pNode1VarTmp->mass_s + pNode1VarTmp->mass_f) * N1Tmp
					+ pNode2VarTmp->dMomentum1_m / (pNode2VarTmp->mass_s + pNode2VarTmp->mass_f) * N2Tmp);
			// update location of particle
			pPclTmp->x += pNode1VarTmp->u1_s * N1Tmp + pNode2VarTmp->u1_s * N2Tmp;
			// update strain increment of particle
			pPclVarTmp->dstrain11 =
				pNode1VarTmp->u1_s * pPclVarTmp->dN1_dx + pNode2VarTmp->u1_s * pPclVarTmp->dN2_dx;
			
			//std::cout << "us: " << pNode1VarTmp->u1_s << " " << pNode2VarTmp->u1_s 
			//	<< " dstrain11: " << pPclVarTmp->dstrain11 << std::endl;
			
			// fluid phase
			pPclTmp->w1 = pNode1VarTmp->w1 * N1Tmp + pNode2VarTmp->w1 * N2Tmp;
			pPclVarTmp->dw_dx =
				pNode1VarTmp->u1_f * pPclVarTmp->dN1_dx + pNode2VarTmp->u1_f * pPclVarTmp->dN2_dx;
			
			//std::cout << "uf: " << pNode1VarTmp->u1_f << " " << pNode2VarTmp->u1_f
			//	<< " dw_dx: " << pPclVarTmp->dw_dx << std::endl;

			//std::cout << "iter_id: " << iteration_index 
			//	<< " dstrain11: " << pPclVarTmp->dstrain11 
			//	<< " dw_dx: " << pPclVarTmp->dw_dx << std::endl;
		}
	}

	return 0;
}


int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	updatePhysicalPropertyAtParticle_SingleObject(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	double dstrain[6], estress[6];
	double destress[6], destrain[6], dpstrain[6];
	double F_det;
	double dp; // increment of pore pressure
	double dDensity_f; // increment of density
	double dstrain_v_f; // volumetric strain of fluid

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->particles + i;
		pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
		if (pPclTmp->isInMesh)
		{
			//std::cout << "pcl_id: " << pPclTmp->index;
			// ------------------- Fluid Phase -------------------
			// Take compression as positive
			dstrain_v_f = -(pPclVarTmp->dw_dx + pPclVarTmp->dstrain11) / pPclTmp->n;
			// update density of fluid phase
			// dV / V = dDensity_f / density_f = dp / Kf;
			dDensity_f = pPclTmp->density_f * dstrain_v_f;
			pPclTmp->density_f += dDensity_f;
			// update pore pressure
			dp = pPclTmp->Kf * dstrain_v_f;
			pPclTmp->p += dp;

			//std::cout << "dstrain_v_f: " << dstrain_v_f 
			//	<< " p: " << pPclTmp->p
			//	<< " dp: " << dp << std::endl;

			// ------------------- solid phase -------------------
			// Assume density of solid phase remain unchanged
			//pPclTmp->density_s = constant;
			// Integrate constitute model
			dstrain[0] = pPclVarTmp->dstrain11;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			estress[0] = pPclTmp->estress11;
			estress[1] = 0.0;
			estress[2] = 0.0;
			estress[3] = 0.0;
			estress[4] = 0.0;
			estress[5] = 0.0;

			pPclTmp->cm->integration_explicit(dstrain, estress, destress, destrain, dpstrain);
			//std::cout << " dstrain: " << pPclVarTmp->dstrain11 << " destress: " << destress[0];
			//std::cout << " p: " << pPclTmp->p << std::endl;

			// Update effective stress
			pPclTmp->estress11 += destress[0];
			// Update total stress
			pPclTmp->stress11 = pPclTmp->estress11 - pPclTmp->p;
			// update strain (also assume that the strain increment is Jaumann rate)
			pPclTmp->strain11 += pPclVarTmp->dstrain11;
			pPclTmp->estrain11 += destrain[0];
			pPclTmp->pstrain11 += dpstrain[0];

			//std::cout << " str: " << pPclTmp->estress11
			//	<< " dstr: " << destress[0] << std::endl;

			//std::cout << "n: " << pPclTmp->n
			//	<< " dstrain_v_f: " << dstrain_v_f
			//	<< " dp: " << dp
			//	<< " dstr: " << destress[0] << std::endl;
			//std::cout << "n: " << pPclTmp->n
			//	<< " dstrain_v_f: " << dstrain_v_f
			//	<< " p: " << pPclTmp->p
			//	<< " str: " << pPclTmp->estress11 << std::endl;

			// update porosity
			F_det = 1.0 + pPclVarTmp->dstrain11;
			pPclTmp->n = 1.0 - (1.0 - pPclTmp->n) / F_det;
			// update mass of fluid phase
			pPclTmp->mass_f = pPclTmp->density_f * pPclTmp->n *
				(pPclTmp->mass_s / pPclTmp->density_s / (1.0 - pPclTmp->n));

			//std::cout << pPclTmp->mass_f / pPclTmp->density_f /(pPclTmp->mass_s / pPclTmp->density_s + pPclTmp->mass_f / pPclTmp->density_f) << std::endl;
		}
	}

	return 0;
}
