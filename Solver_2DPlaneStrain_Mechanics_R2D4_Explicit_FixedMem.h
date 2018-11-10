#ifndef __SOLVER_2DPLANESTRAIN_MECHANICS_R2D4_EXPLICIT_FIXEDMEM_H__
#define __SOLVER_2DPLANESTRAIN_MECHANICS_R2D4_EXPLICIT_FIXEDMEM_H__

#include <vector>

#include "Solver.h"
#include "Mesh_R2D4.h"
#include "Particle_2D_Mechanics.h"
#include "OutputRequest.h"

class Mesh_R2D4;
class ObjectByParticle;

class Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem : public Solver
{
protected:
	// Background Mesh
	Mesh_R2D4 &mesh;
	double *nodeXCoords;
	double *nodeYCoords;
	Node_R2D *nodes;
	size_t nodeNum;
	Element_R2D4 *elements;
	size_t elementNum;

	// Objects described by material points
	std::vector<ObjectByParticle_2D_Mechanics> &pcl_objects;
	size_t object_num;
	ObjectByParticle_2D_Mechanics *curObject;
	size_t curObjectId; // start from 0

	// point to memory of all nodal variables
	NodeVar_2D_Mechanics *nodeVarMem;
	// point to memeory of all particle variables
	ParticleVar_2D_Mechanics *particleVarMem;

public:
	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(
		const double time_step,
		Mesh_R2D4 &mh,
		std::vector<ObjectByParticle_2D_Mechanics> &pcl_objs,
		OutputRequest &out);
	~Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem();

	// Calculate time step size and store it into dt.
	void CFLTimeStep(void);

	int init(void);

	int iteration(void);

	// Initialize each interations
	int initStep_SingleObject(void);

	// Map mass and linear momentum from node to particle.
	int mapPhysicalProperyToNode_SingleObject(void);

	// calculate nodal force caused by stress within element.
	int calInternalForce_SingleObject(void);

	// Calculate nodal force caused by body force and surface force.
	int calExternalForce_SingleObject(void);

	// Integrate momentum
	// Update velocity
	// Apply boundary conditions
	// Do not detect contact
	int updatePhysicalPropertyAtNode_SingleObject(void);

	int calContactForce(void);

	int mapPhysicalPropertyToParticle_SingleObject(void);

	int updatePhysicalPropertyAtParticle_SingleObject(void);
};

#endif