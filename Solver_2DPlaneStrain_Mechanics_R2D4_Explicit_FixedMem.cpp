#ifdef _DEBUG
#include <iostream>
#endif

#include "Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem.h"

#include "GaussIntegrationMacro.h"

// Solver
Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(	
	const double time_step,
	Mesh_R2D4 &mh,
	std::vector<ObjectByParticle_2D_Mechanics> &pcl_objs,
	OutputRequest &out) :
	mesh(mh), pcl_objects(pcl_objs),
	Solver(time_step, out),
	nodeVarMem(nullptr), particleVarMem(nullptr)
{
	nodeXCoords = mesh.nodeXCoords;
	nodeYCoords = mesh.nodeYCoords;
	nodes = mesh.nodes;
	nodeNum = mesh.nodeNum;
	elements = mesh.elements;
	elementNum = mesh.elementNum;

	object_num = pcl_objects.size();
}

Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	~Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem()
{
	if (nodeVarMem) delete[] nodeVarMem;
	if (particleVarMem) delete[] particleVarMem;
}

void Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::CFLTimeStep(void)
{
	//dt = 0.0; // need further work here...
}

int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::init(void)
{
	size_t i, j, k;
	size_t total_pcl_num, pcl_num;
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	Node_R2D *pNodeTmp;
	NodeVar_2D_Mechanics *pNodeVarTmp;

	// Initialize calculation variables at nodes
	nodeVarMem = new NodeVar_2D_Mechanics[nodeNum * object_num];
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
	particleVarMem = new ParticleVar_2D_Mechanics[total_pcl_num];
	if (!particleVarMem) throw std::exception("Fail to allocate particle variables");

	k = 0;
	for (i = 0; i < object_num; i++)
	{
		curObject = &pcl_objects[i];
		pcl_num = curObject->particleNum;
		//std::cout << pcl_num << std::endl;
		for (j = 0; j < pcl_num; j++)
		{
			pPclTmp = curObject->getParticleById(j + 1);
			//((LinearElasticityModel*)pPclTmp->cm)->showParam();
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


int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::iteration(void)
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
		mapPhysicalPropertyToParticle_SingleObject();
		updatePhysicalPropertyAtParticle_SingleObject();
	}

	return 0;
}


int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	initStep_SingleObject(void)
{
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	Element_R2D4 *pElemTmp;
	Node_R2D *pNodeTmp;
	NodeVar_2D_Mechanics *pNodeVarTmp;
	size_t i;
	size_t pcl_num_tmp;

	// for calculating one order derivatives
	double xi, eta;
	double x1, x2, x3, x4, y1, y2, y3, y4;
	double dN1_dxi, dN1_deta, dN2_dxi, dN2_deta;
	double dN3_dxi, dN3_deta, dN4_dxi, dN4_deta;
	double dx_dxi, dx_deta, dy_dxi, dy_deta;
	double dxi_dx, dxi_dy, deta_dx, deta_dy;
	double Jdet;

	// init all calculation variables
	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		pNodeVarTmp->needCal = false;
	}

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);
			//std::cout << pPclVarTmp->particle->index << std::endl;

			pElemTmp = mesh.findInWhichElement(pPclTmp->x, pPclTmp->y,
				static_cast<Element_R2D4 *>(pPclVarTmp->element));
			//std::cout << pPclTmp->index << ": " << pElemTmp->index << std::endl;
			if (pElemTmp) // if particles still lie in the mesh
			{
				pPclVarTmp->element = pElemTmp;

				// ---------------- Calculate natural coordinate -----------------
				mesh.calNaturalCoords(pElemTmp, pPclTmp->x, pPclTmp->y, &xi, &eta);
				pPclVarTmp->xi = xi;
				pPclVarTmp->eta = eta;
				//std::cout << "p_id: " << pPclTmp->index << std::endl;
				//std::cout << "x: " << pPclTmp->x << "; y: " << pPclTmp->y << std::endl;
				//std::cout << "xi: " << pPclVarTmp->xi << "; eta: " << pPclVarTmp->eta << std::endl;

				// ------------------ Calculate shape function ------------------- 
				pPclVarTmp->N1 = Mesh_R2D4::N1(xi, eta);
				pPclVarTmp->N2 = Mesh_R2D4::N2(xi, eta);
				pPclVarTmp->N3 = Mesh_R2D4::N3(xi, eta);
				pPclVarTmp->N4 = Mesh_R2D4::N4(xi, eta);
				//std::cout << "N1: " << pPclVarTmp->N1 << "; N2: " << pPclVarTmp->N2
				//		<< "; N3: " << pPclVarTmp->N3 << "; N4: " << pPclVarTmp->N4 << std::endl;

				// ---- Calculate one order derivatives of the shape function ----
				// calculate Jacobian matrix
				dN1_dxi  = Mesh_R2D4::dN1_dxi(xi, eta);
				dN1_deta = Mesh_R2D4::dN1_deta(xi, eta);
				dN2_dxi  = Mesh_R2D4::dN2_dxi(xi, eta);
				dN2_deta = Mesh_R2D4::dN2_deta(xi, eta);
				dN3_dxi  = Mesh_R2D4::dN3_dxi(xi, eta);
				dN3_deta = Mesh_R2D4::dN3_deta(xi, eta);
				dN4_dxi  = Mesh_R2D4::dN4_dxi(xi, eta);
				dN4_deta = Mesh_R2D4::dN4_deta(xi, eta);

				/* 
				The Jacobian matrix:
					dx_dxi, dx_deta,
					dy_dxi, dy_deta,
				*/
				x1 = x4 = nodeXCoords[pElemTmp->xIndex];
				x2 = x3 = nodeXCoords[pElemTmp->xIndex + 1];
				y1 = y2 = nodeYCoords[pElemTmp->yIndex];
				y3 = y4 = nodeYCoords[pElemTmp->yIndex + 1];
				dx_dxi  = dN1_dxi  * x1 + dN2_dxi  * x2 + dN3_dxi  * x3 + dN4_dxi  * x4;
				dx_deta = dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3 + dN4_deta * x4;
				dy_dxi  = dN1_dxi  * y1 + dN2_dxi  * y2 + dN3_dxi  * y3 + dN4_dxi  * y4;
				dy_deta = dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3 + dN4_deta * y4;
				//std::cout << "dx_dxi: " << dx_dxi << "; dx_deta: " << dx_deta
				//		<< "; dy_dxi: " << dy_dxi << "; dy_deta: " << dy_deta << std::endl;
				
				// determinant of Jacobian matrix
				Jdet = dx_dxi * dy_deta - dx_deta * dy_dxi;
				/* The inverse of Jacobian matrix:
				dxi_dx,  dxi_dy,
				deta_dx, deta_dy, */
				dxi_dx  =  dy_deta / Jdet;
				dxi_dy  = -dx_deta / Jdet;
				deta_dx = -dy_dxi  / Jdet;
				deta_dy =  dx_dxi  / Jdet;
				//std::cout << "dxi_dx: " << dxi_dx << "; dxi_dy: " << dxi_dy
				//	<< "; deta_dx: " << deta_dx << "; deta_dy: " << deta_dy << std::endl;

				pPclVarTmp->dN1_dx = dN1_dxi * dxi_dx + dN1_deta * deta_dx;
				pPclVarTmp->dN1_dy = dN1_dxi * dxi_dy + dN1_deta * deta_dy;
				pPclVarTmp->dN2_dx = dN2_dxi * dxi_dx + dN2_deta * deta_dx;
				pPclVarTmp->dN2_dy = dN2_dxi * dxi_dy + dN2_deta * deta_dy;
				pPclVarTmp->dN3_dx = dN3_dxi * dxi_dx + dN3_deta * deta_dx;
				pPclVarTmp->dN3_dy = dN3_dxi * dxi_dy + dN3_deta * deta_dy;
				pPclVarTmp->dN4_dx = dN4_dxi * dxi_dx + dN4_deta * deta_dx;
				pPclVarTmp->dN4_dy = dN4_dxi * dxi_dy + dN4_deta * deta_dy;
				/*
				std::cout << "dN1_dx: " << pPclVarTmp->dN1_dx
					<< "; dN1_dy: " << pPclVarTmp->dN1_dy
					<< "; dN2_dx: " << pPclVarTmp->dN2_dx
					<< "; dN2_dy: " << pPclVarTmp->dN2_dy << std::endl;
				std::cout << "dN3_dx: " << pPclVarTmp->dN3_dx
					<< "; dN3_dy: " << pPclVarTmp->dN3_dy
					<< "; dN4_dx: " << pPclVarTmp->dN4_dx
					<< "; dN4_dy: " << pPclVarTmp->dN4_dy << std::endl;
				*/

				/* ---------------------------------------------------------------
					Initialize nodal variables
				---------------------------------------------------------------- */
#define			INIT_NODEVAR(pNodeVar) \
				do { \
					if (!(pNodeVar->needCal)) \
					{ \
						pNodeVar->needCal = true; \
						pNodeVar->mass = 0.0; \
						pNodeVar->momentum1 = 0.0; \
						pNodeVar->momentum2 = 0.0; \
						pNodeVar->internalForce1 = 0.0; \
						pNodeVar->internalForce2 = 0.0; \
						pNodeVar->externalForce1 = 0.0; \
						pNodeVar->externalForce2 = 0.0; \
						pNodeVar->contactForce1 = 0.0; \
						pNodeVar->contactForce2 = 0.0; \
					} \
				} while(0)

				// node 1
				pNodeTmp = pElemTmp->node1;
				pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar1 = pNodeVarTmp;
				//std::cout << pNodeTmp->index << std::endl;

				// node 2
				pNodeTmp = pElemTmp->node2;
				pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar2 = pNodeVarTmp;
				//std::cout << pNodeTmp->index << std::endl;

				// node 3
				pNodeTmp = pElemTmp->node3;
				pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar3 = pNodeVarTmp;
				//std::cout << pNodeTmp->index << std::endl;

				// node 4
				pNodeTmp = pElemTmp->node4;
				pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
				INIT_NODEVAR(pNodeVarTmp);
				pPclVarTmp->nodeVar4 = pNodeVarTmp;
				//std::cout << pNodeTmp->index << std::endl;
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
	1. finds in which element each particle lies (and whether they lie outside the mesh);
	2. allocates space for nodeVar on each nodes that are needed for calculation.
*/
int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	mapPhysicalProperyToNode_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	NodeVar_2D_Mechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pPclVarTmp->nodeVar1);
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N1;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N1;
			pNodeVarTmp->momentum2 += pPclTmp->momentum2 * pPclVarTmp->N1;
			// node 2
			pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pPclVarTmp->nodeVar2);
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N2;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N2;
			pNodeVarTmp->momentum2 += pPclTmp->momentum2 * pPclVarTmp->N2;
			// node 3
			pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pPclVarTmp->nodeVar3);
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N3;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N3;
			pNodeVarTmp->momentum2 += pPclTmp->momentum2 * pPclVarTmp->N3;
			// node 4
			pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pPclVarTmp->nodeVar4);
			pNodeVarTmp->mass += pPclTmp->mass * pPclVarTmp->N4;
			pNodeVarTmp->momentum1 += pPclTmp->momentum1 * pPclVarTmp->N4;
			pNodeVarTmp->momentum2 += pPclTmp->momentum2 * pPclVarTmp->N4;
		}
	}

	/*
	for (i = 0; i < mesh.nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		std::cout << "n_id: " << nodes[i].index
			<< "; m: " << pNodeVarTmp->mass << std::endl;
		std::cout << "; m1: " << pNodeVarTmp->momentum1
			<< "; m2: " << pNodeVarTmp->momentum2 << std::endl;
	}
	*/

	return 0;
}


int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	calInternalForce_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	NodeVar_2D_Mechanics *pNodeVarTmp;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			//std::cout << "s11: " << pPclTmp->stress11 << "; s12: " << pPclTmp->stress12
			//		<< "; s22: " << pPclTmp->stress22 << std::endl;
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);
			// node 1
			pNodeVarTmp = pPclVarTmp->nodeVar1;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN1_dx * pPclTmp->stress11 + pPclVarTmp->dN1_dy * pPclTmp->stress12);
			pNodeVarTmp->internalForce2 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN1_dx * pPclTmp->stress12 + pPclVarTmp->dN1_dy * pPclTmp->stress22);
			// node 2
			pNodeVarTmp = pPclVarTmp->nodeVar2;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN2_dx * pPclTmp->stress11 + pPclVarTmp->dN2_dy * pPclTmp->stress12);
			pNodeVarTmp->internalForce2 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN2_dx * pPclTmp->stress12 + pPclVarTmp->dN2_dy * pPclTmp->stress22);
			// node 3
			pNodeVarTmp = pPclVarTmp->nodeVar3;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN3_dx *pPclTmp->stress11 + pPclVarTmp->dN3_dy * pPclTmp->stress12);
			pNodeVarTmp->internalForce2 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN3_dx *pPclTmp->stress12 + pPclVarTmp->dN3_dy * pPclTmp->stress22);
			// node 4
			pNodeVarTmp = pPclVarTmp->nodeVar4;
			pNodeVarTmp->internalForce1 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN4_dx *pPclTmp->stress11 + pPclVarTmp->dN4_dy * pPclTmp->stress12);
			pNodeVarTmp->internalForce2 += pPclTmp->mass / pPclTmp->density *
				(pPclVarTmp->dN4_dx *pPclTmp->stress12 + pPclVarTmp->dN4_dy * pPclTmp->stress22);
		}
	}

	/*
	for (i = 0; i < mesh.nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		std::cout << "n_id: " << nodes[i].index << std::endl;
		std::cout << "int_f1: " << pNodeVarTmp->internalForce1
				<< "; int_f2: " << pNodeVarTmp->internalForce2 << std::endl;
	}
	*/

	return 0;
}

int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	calExternalForce_SingleObject(void)
{
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	double mass_force_tmp;
	MassForceBC *mfp_iter;
	double surface_force_tmp;
	SurfaceForceBC *sfp_iter;

	// ------------------- Boundary conditions specified on particles ------------------- 
	// body force per unit mass specified on particles 
	for (mfp_iter = curObject->massForceBCs_mem.get_first(); mfp_iter;
		 mfp_iter = curObject->massForceBCs_mem.get_next(mfp_iter))
	{
		pPclTmp = curObject->getParticleById(mfp_iter->index);
		//std::cout << "bodyf1: " << mfp_iter->bodyForce1(t_cal) << " bodyf2: " << mfp_iter->bodyForce2(t_cal) << std::endl;
		if (pPclTmp->isInMesh)
		{
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);
			mass_force_tmp = mfp_iter->massForce(t_cal);
			switch (mfp_iter->dof)
			{
			case DegreeOfFreedom::x:
				// node 1
				pPclVarTmp->nodeVar1->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N1;
				// node 2
				pPclVarTmp->nodeVar2->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N2;
				// node 3
				pPclVarTmp->nodeVar3->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N3;
				// node 4
				pPclVarTmp->nodeVar4->externalForce1 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N4;
				break;
			case DegreeOfFreedom::y:
				// node 1
				pPclVarTmp->nodeVar1->externalForce2 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N1;
				// node 2
				pPclVarTmp->nodeVar2->externalForce2 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N2;
				// node 3
				pPclVarTmp->nodeVar3->externalForce2 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N3;
				// node 4
				pPclVarTmp->nodeVar4->externalForce2 += pPclTmp->mass * mass_force_tmp * pPclVarTmp->N4;
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
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);
			surface_force_tmp = sfp_iter->surfaceForce(t_cal);
			switch (sfp_iter->dof)
			{
			case DegreeOfFreedom::x:
				// node 1
				pPclVarTmp->nodeVar1->externalForce1 += surface_force_tmp * pPclVarTmp->N1;
				// node 2
				pPclVarTmp->nodeVar2->externalForce1 += surface_force_tmp * pPclVarTmp->N2;
				// node 3
				pPclVarTmp->nodeVar3->externalForce1 += surface_force_tmp * pPclVarTmp->N3;
				// node 4
				pPclVarTmp->nodeVar4->externalForce1 += surface_force_tmp * pPclVarTmp->N4;
				break;
			case DegreeOfFreedom::y:
				// node 1
				pPclVarTmp->nodeVar1->externalForce2 += surface_force_tmp * pPclVarTmp->N1;
				// node 2
				pPclVarTmp->nodeVar2->externalForce2 += surface_force_tmp * pPclVarTmp->N2;
				// node 3
				pPclVarTmp->nodeVar3->externalForce2 += surface_force_tmp * pPclVarTmp->N3;
				// node 4
				pPclVarTmp->nodeVar4->externalForce2 += surface_force_tmp * pPclVarTmp->N4;
				break;
			default:
				break;
			}
		}
	}

	return 0;
}

int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	updatePhysicalPropertyAtNode_SingleObject(void)
{
	size_t i;
	Node_R2D *pNodeTmp;
	NodeVar_2D_Mechanics *pNodeVarTmp;
	AccelerationBC *a_iter;
	VelocityBC *v_iter;

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			// calculate nodal force
			pNodeVarTmp->nodalForce1 = pNodeVarTmp->externalForce1
				- pNodeVarTmp->internalForce1 + pNodeVarTmp->contactForce1;
			pNodeVarTmp->nodalForce2 = pNodeVarTmp->externalForce2
				- pNodeVarTmp->internalForce2 + pNodeVarTmp->contactForce2;
			// calculate increment of linear momentum
			pNodeVarTmp->dMomentum1 = pNodeVarTmp->nodalForce1 * t_increment_a;
			pNodeVarTmp->dMomentum2 = pNodeVarTmp->nodalForce2 * t_increment_a;
			
			//std::cout << "n_id: " << nodes[i].index << std::endl;
			//std::cout << "dm1: " << pNodeVarTmp->dMomentum1
			//		<< "; dm2: " << pNodeVarTmp->dMomentum2 << std::endl;
			/*if (i == 37 && t_cal > 8.0 && t_cal < 8.1)
			{
				std::cout << "n_id: " << nodes[i].index;
				std::cout << " time: " << t_cal;
				//std::cout << " dmv1: " << pNodeVarTmp->dMomentum1 << " dmv2: " << pNodeVarTmp->dMomentum2 << std::endl;
				std::cout << " ef1: " << pNodeVarTmp->externalForce1 << " ef2: " << pNodeVarTmp->externalForce2;
				std::cout << " if1: " << pNodeVarTmp->internalForce1 << " if2: " << pNodeVarTmp->internalForce2;
				std::cout << " cf1: " << pNodeVarTmp->contactForce1 << " cf2: " << pNodeVarTmp->contactForce2 << std::endl;
				//std::cout << " nf1: " << pNodeVarTmp->nodalForce1 << " nf2: " << pNodeVarTmp->nodalForce2 << std::endl;
			}*/
		}
	}

	// Apply acceleration boundary conditions
	for (a_iter = mesh.accelerationBCs_mem.get_first(); a_iter;
		 a_iter = mesh.accelerationBCs_mem.get_next(a_iter))
	{
		pNodeTmp = mesh.getNodeById(a_iter->index);
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			switch (a_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNodeVarTmp->dMomentum1 = a_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
				break;
			case DegreeOfFreedom::y:
				pNodeVarTmp->dMomentum2 = a_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
				break;
			default:
				break;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			// Integrate linear momentum
			pNodeVarTmp->momentum1 += pNodeVarTmp->dMomentum1;
			pNodeVarTmp->momentum2 += pNodeVarTmp->dMomentum2;
		
			/*if (i == 37 && t_cal > 8.0 && t_cal < 8.1)
			{
				std::cout << "n_id: " << nodes[i].index;
				std::cout << " time: " << t_cal;
				std::cout << " dmv1: " << pNodeVarTmp->dMomentum1 << " dmv2: " << pNodeVarTmp->dMomentum2 << std::endl;
			}*/
		}
	}

	// Apply velocity boundary conditions
	for (v_iter = mesh.velocityBCs_mem.get_first(); v_iter;
		 v_iter = mesh.velocityBCs_mem.get_next(v_iter))
	{
		pNodeTmp = mesh.getNodeById(v_iter->index);
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(pNodeTmp->nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			switch (v_iter->dof)
			{
			case DegreeOfFreedom::x:
				pNodeVarTmp->dMomentum1 = v_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
				pNodeVarTmp->momentum1 = v_iter->v(t_cal) * pNodeVarTmp->mass;
				break;
			case DegreeOfFreedom::y:
				pNodeVarTmp->dMomentum2 = v_iter->a(t_cal) * pNodeVarTmp->mass * t_increment_a;
				pNodeVarTmp->momentum2 = v_iter->v(t_cal) * pNodeVarTmp->mass;
				break;
			default:
				break;
			}
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		pNodeVarTmp = static_cast<NodeVar_2D_Mechanics *>(nodes[i].nodeVar) + curObjectId;
		if (pNodeVarTmp->needCal)
		{
			// calculate velocity
			pNodeVarTmp->v1 = pNodeVarTmp->momentum1 / pNodeVarTmp->mass;
			pNodeVarTmp->v2 = pNodeVarTmp->momentum2 / pNodeVarTmp->mass;
			// calculate increment of displacement
			pNodeVarTmp->u1 = pNodeVarTmp->v1 * t_increment;
			pNodeVarTmp->u2 = pNodeVarTmp->v2 * t_increment;
			
			/*if (i == 37)
			{
				//<< "n_id: " << nodes[i].index 
				std::cout << " time: " << t_cal;
				std::cout << " v1: " << pNodeVarTmp->v1 << " v2: " << pNodeVarTmp->v2;
				std::cout << " mv1: " << pNodeVarTmp->momentum1 << " m: " << pNodeVarTmp->mass << std::endl;
			}*/
		}
	}
	
	//std::cout << "Time: " << t_cal << " mv: " << 

	return 0;
}

int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	calContactForce(void)
{
	return 0;
}

// map momentum, displacement and strain
int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	mapPhysicalPropertyToParticle_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	NodeVar_2D_Mechanics *pNode1VarTmp, *pNode2VarTmp;
	NodeVar_2D_Mechanics *pNode3VarTmp, *pNode4VarTmp;
	double N1Tmp, N2Tmp, N3Tmp, N4Tmp;
	
	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);
		if (pPclTmp->isInMesh)
		{
			//std::cout << "p_id: " << pPclTmp->index << std::endl;
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);

			// --------------------------- node 1 ------------------------------
			pNode1VarTmp = pPclVarTmp->nodeVar1;
			N1Tmp = pPclVarTmp->N1;
			// --------------------------- node 2 ------------------------------
			pNode2VarTmp = pPclVarTmp->nodeVar2;
			N2Tmp = pPclVarTmp->N2;
			// --------------------------- node 3 ------------------------------
			pNode3VarTmp = pPclVarTmp->nodeVar3;
			N3Tmp = pPclVarTmp->N3;
			// --------------------------- node 4 ------------------------------
			pNode4VarTmp = pPclVarTmp->nodeVar4;
			N4Tmp = pPclVarTmp->N4;
			//std::cout << "N1: " << N1Tmp << "; N2: " << N2Tmp << "; N3: " << N3Tmp
			//	<< "; N4: " << N4Tmp << std::endl;

			// update linear momentum of particle
			pPclTmp->momentum1 += pPclTmp->mass *
				(pNode1VarTmp->dMomentum1 / pNode1VarTmp->mass * N1Tmp +
				pNode2VarTmp->dMomentum1 / pNode2VarTmp->mass * N2Tmp +
				pNode3VarTmp->dMomentum1 / pNode3VarTmp->mass * N3Tmp +
				pNode4VarTmp->dMomentum1 / pNode4VarTmp->mass * N4Tmp);
			pPclTmp->momentum2 += pPclTmp->mass *
				(pNode1VarTmp->dMomentum2 / pNode1VarTmp->mass * N1Tmp +
				pNode2VarTmp->dMomentum2 / pNode2VarTmp->mass * N2Tmp +
				pNode3VarTmp->dMomentum2 / pNode3VarTmp->mass * N3Tmp +
				pNode4VarTmp->dMomentum2 / pNode4VarTmp->mass * N4Tmp);

			// update location of particle
			pPclTmp->x += pNode1VarTmp->u1 * N1Tmp + pNode2VarTmp->u1 * N2Tmp
						+ pNode3VarTmp->u1 * N3Tmp + pNode4VarTmp->u1 * N4Tmp;
			pPclTmp->y += pNode1VarTmp->u2 * N1Tmp + pNode2VarTmp->u2 * N2Tmp
						+ pNode3VarTmp->u2 * N3Tmp + pNode4VarTmp->u2 * N4Tmp;

			/*if (i == 6 && t_cal > 8.0 && t_cal < 8.1)
			{
				std::cout << "Time: " << t_cal << " n1: " << pNode1VarTmp->u2 << " n2: " << pNode2VarTmp->u2
					<< " n3: " << pNode3VarTmp->u2 << " n4: " << pNode4VarTmp->u2 << std::endl;
				std::cout << pNode3VarTmp->node->index << std::endl;
			}*/

			// update strain increment of particle
			pPclVarTmp->dstrain11 =
				  pNode1VarTmp->u1 * pPclVarTmp->dN1_dx
				+ pNode2VarTmp->u1 * pPclVarTmp->dN2_dx
				+ pNode3VarTmp->u1 * pPclVarTmp->dN3_dx
				+ pNode4VarTmp->u1 * pPclVarTmp->dN4_dx;
			pPclVarTmp->dstrain22 =
				  pNode1VarTmp->u2 * pPclVarTmp->dN1_dy
				+ pNode2VarTmp->u2 * pPclVarTmp->dN2_dy
				+ pNode3VarTmp->u2 * pPclVarTmp->dN3_dy
				+ pNode4VarTmp->u2 * pPclVarTmp->dN4_dy;
			pPclVarTmp->dstrain12 = 0.5 *
				 (pNode1VarTmp->u1 * pPclVarTmp->dN1_dy
				+ pNode2VarTmp->u1 * pPclVarTmp->dN2_dy
				+ pNode3VarTmp->u1 * pPclVarTmp->dN3_dy
				+ pNode4VarTmp->u1 * pPclVarTmp->dN4_dy
				+ pNode1VarTmp->u2 * pPclVarTmp->dN1_dx
				+ pNode2VarTmp->u2 * pPclVarTmp->dN2_dx
				+ pNode3VarTmp->u2 * pPclVarTmp->dN3_dx
				+ pNode4VarTmp->u2 * pPclVarTmp->dN4_dx);

			//std::cout << "strain11: " << pPclVarTmp->dstrain11
			//		  << " strain22: " << pPclVarTmp->dstrain22
			//		  << " strain12: " << pPclVarTmp->dstrain12;

			// calculate W for Jaumann rate
			pPclVarTmp->dW12 = 0.5 *
				 (pNode1VarTmp->u1 * pPclVarTmp->dN1_dy
				+ pNode2VarTmp->u1 * pPclVarTmp->dN2_dy
				+ pNode3VarTmp->u1 * pPclVarTmp->dN3_dy
				+ pNode4VarTmp->u1 * pPclVarTmp->dN4_dy
				- pNode1VarTmp->u2 * pPclVarTmp->dN1_dx
				- pNode2VarTmp->u2 * pPclVarTmp->dN2_dx
				- pNode3VarTmp->u2 * pPclVarTmp->dN3_dx
				- pNode4VarTmp->u2 * pPclVarTmp->dN4_dx);
			//std::cout << "; dW12: " << pPclVarTmp->dW12 << std::endl;
		}
	}

	return 0;
}


int Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem::
	updatePhysicalPropertyAtParticle_SingleObject(void)
{
	size_t i;
	size_t pcl_num_tmp;
	Particle_2D_Mechanics *pPclTmp;
	ParticleVar_2D_Mechanics *pPclVarTmp;
	double dstrain[6], destrain[6], dpstrain[6];
	double dstress[6];
	double F_det;

	pcl_num_tmp = curObject->particleNum;
	for (i = 0; i < pcl_num_tmp; i++)
	{
		pPclTmp = curObject->getParticleById(i + 1);

		if (pPclTmp->isInMesh) // only update element when in mesh
		{
			pPclVarTmp = static_cast<ParticleVar_2D_Mechanics *>(pPclTmp->particleVar);

			// update density
			F_det = (1.0 + pPclVarTmp->dstrain11 + pPclVarTmp->dstrain22);
			pPclTmp->density /= F_det;

			// integrate constitute
			dstrain[0] = pPclVarTmp->dstrain11;
			dstrain[1] = pPclVarTmp->dstrain22;
			dstrain[2] = 0.0; // plane strain condition
			dstrain[3] = pPclVarTmp->dstrain12;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;

			//std::cout << "p_id: " << pPclTmp->index << std::endl;
			//((LinearElasticityModel*)pPclTmp->cm)->showParam();
			pPclTmp->cm->integration_explicit(dstrain, pPclTmp->stress, dstress, destrain, dpstrain);

			//std::cout << "strain11: " << pPclVarTmp->dstrain11
			//	<< "; strain12: " << pPclVarTmp->dstrain12
			//	<< "; strain22: " << pPclVarTmp->dstrain22 << std::endl;
			//std::cout << "s11: " << dstress[0] << "; s12: " << dstress[3]
			//		<< "; s22: " << dstress[1] << std::endl;

			/*
			Rotate as Jaumann rate:
				tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
			*/
			dstress[0] +=  pPclVarTmp->dW12 * pPclTmp->stress12 * 2.0;
			dstress[1] += -pPclVarTmp->dW12 * pPclTmp->stress12 * 2.0;
			dstress[3] +=  pPclVarTmp->dW12 * (pPclTmp->stress22 - pPclTmp->stress11);

			//std::cout << "s11: " << dstress[0] << "; s12: " << dstress[3]
			//		<< "; s22: " << dstress[1] << std::endl;

			//update stress
			pPclTmp->stress11 += dstress[0];
			pPclTmp->stress22 += dstress[1];
			pPclTmp->stress33 += dstress[2];
			pPclTmp->stress12 += dstress[3];
			pPclTmp->stress23 += dstress[4];
			pPclTmp->stress31 += dstress[5];

			/*if (i == 6 && t_cal > 8.0 && t_cal < 8.1)
			{
				std::cout << "time: " << t_cal
					<< " dsa11: " << pPclVarTmp->dstrain11
					<< " dsa22: " << pPclVarTmp->dstrain22
					<< " dsa12: " << pPclVarTmp->dstrain12 << std::endl;
			}*/

			// update strain (also assume that strain increment is Jaumann rate)
			pPclVarTmp->dstrain11 +=  pPclVarTmp->dW12 * pPclTmp->strain12 * 2.0;
			pPclVarTmp->dstrain22 += -pPclVarTmp->dW12 * pPclTmp->strain12 * 2.0;
			pPclVarTmp->dstrain12 +=  pPclVarTmp->dW12 * (pPclTmp->strain22 - pPclTmp->strain11);
			pPclTmp->strain11 += pPclVarTmp->dstrain11;
			pPclTmp->strain22 += pPclVarTmp->dstrain22;
			pPclTmp->strain12 += pPclVarTmp->dstrain12;

			/*if (i == 6 && t_cal > 8.0 && t_cal < 8.1)
			{
				std::cout << "time: " << t_cal 
					<< " dsa11: " << pPclVarTmp->dstrain11 
					<< " dsa22: " << pPclVarTmp->dstrain22 
					<< " dsa12: " << pPclVarTmp->dstrain12 << std::endl;
			}*/

			destrain[0] +=  pPclVarTmp->dW12 * pPclTmp->estrain12 * 2.0;
			destrain[1] += -pPclVarTmp->dW12 * pPclTmp->estrain12 * 2.0;
			destrain[3] +=  pPclVarTmp->dW12 * (pPclTmp->estrain22 - pPclTmp->estrain11);
			pPclTmp->estrain11 += destrain[0];
			pPclTmp->estrain22 += destrain[1];
			pPclTmp->estrain12 += destrain[3];

			dpstrain[0] +=  pPclVarTmp->dW12 * pPclTmp->pstrain12 * 2.0;
			dpstrain[1] += -pPclVarTmp->dW12 * pPclTmp->pstrain12 * 2.0;
			dpstrain[3] +=  pPclVarTmp->dW12 * (pPclTmp->pstrain22 - pPclTmp->pstrain11);
			pPclTmp->pstrain11 += dpstrain[0];
			pPclTmp->pstrain22 += dpstrain[1];
			pPclTmp->pstrain12 += dpstrain[3];
		}
	}

	return 0;
}
