#ifndef __SOLVER_1D_MECHANICS_1D2_EXPLICIT_H__
#define __SOLVER_1D_MECHANICS_1D2_EXPLICIT_H__

#include <vector>
#include "MemoryPoolTemplate.hpp"
#include "Data_1D_Mechanics_1D2.h"
#include "OutputRequest.h"

class Solver_1D_Mechanics_1D2_Explicit : public Solver
{
protected:
	// Data: contain all model data
	Data_1D_Mechanics_1D2 &data;

	// Background Mesh
	Mesh_1D2 &mesh;
	double *nodeCoords;
	Node_1D *nodes;
	size_t nodeNum;
	Element_1D2 *elements;
	size_t elementNum;

	// Objects described by material points
	std::vector<ObjectByParticle_1D_Mechanics> &objects;
	size_t object_num;
	ObjectByParticle_1D_Mechanics *curObject;

	// Memory pools of NodeVar and ParticleVar
	MemoryPoolTemplate<ParticleVar_1D_Mechanics> particleVarPool;
	MemoryPoolTemplate<NodeVar_1D_Mechanics> nodeVarPool;

	// Index of node where contact may be possible
	std::vector<size_t> node_contact;

public:
	Solver_1D_Mechanics_1D2_Explicit(Data_1D_Mechanics_1D2 &da, OutputRequest &out);
	~Solver_1D_Mechanics_1D2_Explicit();
	
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
	int updatePhysicalPropertyAtNode_noContact(void);

	int calContactForce(void);

	int mapPhysicalPropertyToParticle_SingleObject(void);

	int updatePhysicalPropertyAtParticle_SingleObject(void);
};

#endif