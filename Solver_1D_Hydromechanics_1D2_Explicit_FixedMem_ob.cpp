#ifdef _DEBUG
#include <iostream>
#endif

#include "Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_ob.h"

#include "GaussIntegrationMacro.h"

Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(
	const double time_step,
	Mesh_1D2 &mh,
	std::vector<ObjectByParticle_1D_Hydromechanics> &pcl_objs,
	OutputRequest &out) :
	mesh(mh), pcl_objects(pcl_objs),
	Solver(time_step, &out),
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

	// Initialize calculation variables at nodes
	nodeVarMem = new NodeVar_1D_Hydromechanics[nodeNum * object_num];
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
	particleVarMem = new ParticleVar_1D_Hydromechanics[total_pcl_num];
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
			pPclVarTmp->v1_s = 0.0;
			pPclVarTmp->v1_f = 0.0;
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
		calInternalAndSeepageForce_SingleObject();
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
int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::initStep_SingleObject(void)
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
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
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
				
				/*std::cout << "pcl_id: " << pPclTmp->index 
					<< " dN1_dxi: "<< pPclVarTmp->dN1_dx
					<< " dN2_dxi: " << pPclVarTmp->dN2_dx
					<< " n: " << pPclTmp->n <<std::endl;*/

				pPclVarTmp->unit_weight_fluid = 9.8 * pPclTmp->density_f;
				pPclVarTmp->avgdensity_s = (1.0 - pPclTmp->n) * pPclTmp->density_s;
				pPclVarTmp->avgdensity_f = pPclTmp->n * pPclTmp->density_f;

				/* ---------------------------------------------------------------
				Initialize nodal variables
				---------------------------------------------------------------- */
#define			INIT_NODEVAR(pNodeVar) \
				do { \
					if (!(pNodeVar->needCal)) \
					{ \
						pNodeVarTmp->needCal = true; \
						pNodeVarTmp->mass_s = 0.0; \
						pNodeVarTmp->momentum1_s = 0.0; \
						pNodeVarTmp->internalForce1_s = 0.0; \
						pNodeVarTmp->externalForce1_s = 0.0; \
						pNodeVarTmp->contactForce1_s = 0.0; \
						pNodeVarTmp->mass_f = 0.0; \
						pNodeVarTmp->momentum1_f = 0.0; \
						pNodeVarTmp->internalForce1_f = 0.0; \
						pNodeVarTmp->externalForce1_f = 0.0; \
						pNodeVarTmp->contactForce1_f = 0.0; \
						pNodeVarTmp->seepageForce1 = 0.0; \
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
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			pNodeVarTmp->mass_s += pPclTmp->mass_s * pPclVarTmp->N1;
			pNodeVarTmp->momentum1_s += pPclTmp->momentum1_s * pPclVarTmp->N1;
			pNodeVarTmp->mass_f += pPclTmp->mass_f * pPclVarTmp->N1;
			pNodeVarTmp->momentum1_f += pPclTmp->momentum1_f * pPclVarTmp->N1;
			// node 2
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			pNodeVarTmp->mass_s += pPclTmp->mass_s * pPclVarTmp->N2;
			pNodeVarTmp->momentum1_s += pPclTmp->momentum1_s * pPclVarTmp->N2;
			pNodeVarTmp->mass_f += pPclTmp->mass_f * pPclVarTmp->N2;
			pNodeVarTmp->momentum1_f += pPclTmp->momentum1_f * pPclVarTmp->N2;
		}
	}

	/*
	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_1D_Hydromechanics *>(nodes[i].nodeVar);
		std::cout << "n_id: " << pNodeVarTmp->node->index
			<< " m_s: " << pNodeVarTmp->mass_s << " mv_s: " << pNodeVarTmp->momentum1_s
			<< " m_f: " << pNodeVarTmp->mass_f << " mv_f: " << pNodeVarTmp->momentum1_f << std::endl;
	}*/

	return 0;
}


int Solver_1D_Hydromechanics_1D2_Explicit_FixedMem::
	calInternalAndSeepageForce_SingleObject(void)
{
	size_t i, pcl_num_tmp;
	Particle_1D_Hydromechanics *pPclTmp;
	ParticleVar_1D_Hydromechanics *pPclVarTmp;
	NodeVar_1D_Hydromechanics *pNodeVarTmp;
	double n_tmp;
	double volume_s_tmp, volume_f_tmp;
	double seepage_force;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			//std::cout << "p_id: " << pPclTmp->index;

			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			n_tmp = pPclTmp->n;
			volume_s_tmp = pPclTmp->mass_s / pPclVarTmp->avgdensity_s;
			volume_f_tmp = pPclTmp->mass_f / pPclVarTmp->avgdensity_f;
			
			// ------------------- Node 1 ----------------------
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			// solid phase
			pNodeVarTmp->internalForce1_s += pPclVarTmp->dN1_dx * 
				(pPclTmp->estress11 - (1.0 - n_tmp) * pPclTmp->p) * volume_s_tmp;
			// fluid phase
			pNodeVarTmp->internalForce1_f += pPclVarTmp->dN1_dx * 
				(n_tmp * -pPclTmp->p) * volume_f_tmp;
			// seepage force 
			seepage_force = pPclVarTmp->N1 * n_tmp * n_tmp * pPclVarTmp->unit_weight_fluid / pPclTmp->kappa;
			pNodeVarTmp->seepageForce1 += seepage_force * (pPclVarTmp->v1_f * volume_f_tmp - pPclVarTmp->v1_s * volume_s_tmp);
			//std::cout << " sftmp: " << seepage_force;
			//std::cout << " v1_f: " << pPclVarTmp->v1_f << " vf: " << volume_f_tmp;
			//std::cout << " v1_s: " << pPclVarTmp->v1_s << " vs: " << volume_s_tmp;
			//std::cout << " sf: " << pNodeVarTmp->seepageForce1;

			// ------------------- Node 2 ----------------------
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			// solid phase
			pNodeVarTmp->internalForce1_s += pPclVarTmp->dN2_dx *
				(pPclTmp->estress11 - (1.0 - n_tmp) * pPclTmp->p) * volume_s_tmp;
			// fluid phase
			pNodeVarTmp->internalForce1_f += pPclVarTmp->dN2_dx *
				(n_tmp * -pPclTmp->p) * volume_f_tmp;
			// seepage force
			seepage_force = pPclVarTmp->N2 * n_tmp * n_tmp * pPclVarTmp->unit_weight_fluid / pPclTmp->kappa;
			pNodeVarTmp->seepageForce1 += seepage_force * (pPclVarTmp->v1_f * volume_f_tmp - pPclVarTmp->v1_s * volume_s_tmp);
			//std::cout << " sftmp: " << seepage_force;
			//std::cout << " sf: " << pNodeVarTmp->seepageForce1 << std::endl;
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
				// node 1
				pNode1VarTmp->externalForce1_s += pPclVarTmp->N1 * mass_force_tmp * pPclTmp->mass_s;
				// node 2
				pNode2VarTmp->externalForce1_s += pPclVarTmp->N2 * mass_force_tmp * pPclTmp->mass_s;
				break;
			case DegreeOfFreedom::x_f:
				// node 1
				pNode1VarTmp->externalForce1_f += pPclVarTmp->N1 * mass_force_tmp * pPclTmp->mass_f;
				// node 2
				pNode2VarTmp->externalForce1_f += pPclVarTmp->N2 * mass_force_tmp * pPclTmp->mass_f;
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
			switch (sfp_iter->dof) // not sure ?????
			{
			case DegreeOfFreedom::x:
				pNode1VarTmp->externalForce1_s += pPclVarTmp->N1 * surface_force_tmp;
				pNode2VarTmp->externalForce1_s += pPclVarTmp->N2 * surface_force_tmp;
				break;
			case DegreeOfFreedom::x_f:
				pNode1VarTmp->externalForce1_s += -pPclVarTmp->N1 * (1.0 - pPclTmp->n) * surface_force_tmp;
				pNode2VarTmp->externalForce1_s += -pPclVarTmp->N2 * (1.0 - pPclTmp->n) * surface_force_tmp;
				pNode1VarTmp->externalForce1_f += -pPclVarTmp->N1 * pPclTmp->n * surface_force_tmp;
				pNode2VarTmp->externalForce1_f += -pPclVarTmp->N2 * pPclTmp->n * surface_force_tmp;
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
			pNodeVarTmp->nodalForce1_s = pNodeVarTmp->externalForce1_s
				- pNodeVarTmp->internalForce1_s + pNodeVarTmp->seepageForce1
				+ pNodeVarTmp->contactForce1_s;
			pNodeVarTmp->nodalForce1_f = pNodeVarTmp->externalForce1_f
				- pNodeVarTmp->internalForce1_f - pNodeVarTmp->seepageForce1
				+ pNodeVarTmp->contactForce1_f;
			// calculate increment of linear momentum
			pNodeVarTmp->dMomentum1_s = pNodeVarTmp->nodalForce1_s * t_increment_a;
			pNodeVarTmp->dMomentum1_f = pNodeVarTmp->nodalForce1_f * t_increment_a;
		
			/*
			std::cout << "n_id: " << pNodeVarTmp->node->index
				<< " nf_s: " << pNodeVarTmp->nodalForce1_s
				<< " nf_f: " << pNodeVarTmp->nodalForce1_f;

			std::cout << " if_s: " << pNodeVarTmp->internalForce1_s
				<< " if_f: " << pNodeVarTmp->internalForce1_f;

			std::cout << " sf: " << pNodeVarTmp->seepageForce1 << std::endl;
			*/
		
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
				pNodeVarTmp->dMomentum1_s = a_iter->a(t_cal) * pNodeVarTmp->mass_s * t_increment_a;
				break;
			// currently no acceleration boudnary for fluid phase
			case DegreeOfFreedom::x_f:
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
			pNodeVarTmp->momentum1_s += pNodeVarTmp->dMomentum1_s;
			pNodeVarTmp->momentum1_f += pNodeVarTmp->dMomentum1_f;
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
				pNodeVarTmp->dMomentum1_s = v_iter->a(t_cal) * pNodeVarTmp->mass_s * t_increment_a;;
				pNodeVarTmp->momentum1_s = v_iter->v(t_cal) * pNodeVarTmp->mass_s;
				break;
			case DegreeOfFreedom::x_f:
				pNodeVarTmp->dMomentum1_f = v_iter->a(t_cal) * pNodeVarTmp->mass_f * t_increment_a;;
				pNodeVarTmp->momentum1_f = v_iter->v(t_cal) * pNodeVarTmp->mass_f;
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
			pNodeVarTmp->v1_s = pNodeVarTmp->momentum1_s / pNodeVarTmp->mass_s;
			pNodeVarTmp->v1_f = pNodeVarTmp->momentum1_f / pNodeVarTmp->mass_f;
			// calculate increment of displacement
			pNodeVarTmp->u_s = pNodeVarTmp->v1_s * t_increment;
			pNodeVarTmp->u_f = pNodeVarTmp->v1_f * t_increment;
		}
	}

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
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
			N1Tmp = pPclVarTmp->N1;
			N2Tmp = pPclVarTmp->N2;
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			pNode2VarTmp = pPclVarTmp->nodeVar2;

			// update linear momentum of particle
			pPclTmp->momentum1_s += pPclTmp->mass_s *
				 (pNode1VarTmp->v1_s * N1Tmp + pNode2VarTmp->v1_s * N2Tmp);
			// update location of particle
			pPclTmp->x += pNode1VarTmp->u_s * N1Tmp + pNode2VarTmp->u_s * N2Tmp;
			// update strain increment of particle
			pPclVarTmp->dstrain11 =
				pNode1VarTmp->u_s * pPclVarTmp->dN1_dx + pNode2VarTmp->u_s * pPclVarTmp->dN2_dx;
			
			// fluid phase, need to be modified later, xi, N1, N2 needed to be recalculated ...
			pPclTmp->momentum1_f += pPclTmp->mass_f *
				(pNode1VarTmp->v1_f * N1Tmp + pNode2VarTmp->v1_f * N2Tmp);
			pPclVarTmp->dU_dx =
				pNode1VarTmp->u_f * pPclVarTmp->dN1_dx + pNode2VarTmp->u_f * pPclVarTmp->dN2_dx;
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

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		pPclVarTmp = static_cast<ParticleVar_1D_Hydromechanics *>(pPclTmp->particleVar);
		if (pPclTmp->isInMesh)
		{
			//std::cout << "pcl_id: " << pPclTmp->index;
			// ------------------- Fluid Phase -------------------
			// update pore pressure
			// increment of pore pressure
			// dp + Kf / (1 - n) * ((1 - n) * ui,i + n * Ui,i) = 0
			dp = - pPclTmp->Kf / pPclTmp->n
				* ((1.0 - pPclTmp->n) * pPclTmp->strain11 + pPclTmp->n * pPclVarTmp->dU_dx);
			pPclTmp->p += dp;
			//std::cout << " dU_dx: " << pPclVarTmp->dU_dx;
			// update density of fluid phase
			// dDensity_f / density_f = dp / Kf;
			dDensity_f = dp / pPclTmp->Kf * pPclTmp->density_f;
			pPclTmp->density_f += dDensity_f;
			//std::cout << " density_f: " << pPclTmp->density_f;

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

			// update porosity
			F_det = 1.0 + pPclVarTmp->dstrain11;
			pPclTmp->n = 1.0 - (1.0 - pPclTmp->n) / F_det;
			//std::cout << " n: " << pPclTmp->n;
		}
	}

	return 0;
}
