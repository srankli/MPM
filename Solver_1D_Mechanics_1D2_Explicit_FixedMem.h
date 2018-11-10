#ifndef __SOLVER_1D_MECHANICS_1D2_EXPLICIT_FIXEDMEM_H__
#define __SOLVER_1D_MECHANICS_1D2_EXPLICIT_FIXEDMEM_H__

#include <vector>

#include "Solver.h"
#include "Mesh_1D2.h"
#include "Particle_1D_Mechanics.h"
#include "OutputRequest.h"

class Solver_1D_Mechanics_1D2_Explicit_FixedMem : public Solver
{
protected:
	// Background Mesh
	Mesh_1D2 &mesh;
	double *nodeCoords;
	Node_1D *nodes;
	size_t nodeNum;
	Element_1D2 *elements;
	size_t elementNum;

	// Objects described by material points
	std::vector<ObjectByParticle_1D_Mechanics> &pcl_objects;
	size_t object_num;
	ObjectByParticle_1D_Mechanics *curObject;

	// point to memory of all nodal variables
	NodeVar_1D_Mechanics *nodeVarMem;
	// point to memeory of all particle variables
	ParticleVar_1D_Mechanics *particleVarMem;

public:
	Solver_1D_Mechanics_1D2_Explicit_FixedMem(
		const double time_step,
		Mesh_1D2 &mh,
		std::vector<ObjectByParticle_1D_Mechanics> &pcl_objs,
		OutputRequest &out);
	~Solver_1D_Mechanics_1D2_Explicit_FixedMem();

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