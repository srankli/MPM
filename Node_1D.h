#ifndef __NODE_1D_H__
#define __NODE_1D_H__

#include "Mesh.h"

struct NodeVar;

struct Node_1D : public Node
{
	size_t xIndex;
	
	// velocity of "center of mesh"
	double v1_cm;
};

// variables for calculation
struct NodeVar_1D_Mechanics : NodeVar
{
public:
	// Physical properties
	double mass;
	// linear momentum
	double momentum1;
	// increment of linear momentum
	double dMomentum1;
	// velocity
	double v1;
	// displacement
	double u;
	// internal force
	double internalForce1;
	// external force
	double externalForce1;
	//nodal force
	double nodalForce1;
	
	// contact force
	double contactForce1;
	// normal of surface, for contact force calculation
	// = sum(mp * NIp,i)
	double normal1;
};

struct NodeVar_1D_Hydromechanics : NodeVar
{
public:
	// Physical properties
	// ----- solid phase -----
	double mass_s;
	double momentum1_s; // linear momentum
	// increment of linear momentum
	double dMomentum1_s;
	double v1_s; // velocity
	double u_s; // displacement
	
	double internalForce1_s; // internal force
	double externalForce1_s; // external force
	double contactForce1_s; // contact force
	// = externalForce1_s - internalForce1_s + seepageForce1 (+ contactForce_s)
	double nodalForce1_s;

	// ----- fluid phase -----
	double mass_f;
	double momentum1_f;
	double dMomentum1_f;
	double v1_f;
	double u_f;

	double internalForce1_f; // internal force
	double externalForce1_f; // external force
	// = externalForce1_s - internalForce1_s - seepageForce1 (+ contactForce_f)
	double contactForce1_f;
	double nodalForce1_f;

	double seepageForce1; // seepage force
};

#endif