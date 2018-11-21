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
	// fluid - solid mixture
	double momentum1_m;
	// increment of linear momentum
	double dMomentum1_m;
	double externalForce1_m; // external force
	double internalForce1_m; // internal force
	double contactForce1_m;  // contact force
	// = externalForce1_m - internalForce1_m (+ contactForce_m)
	double nodalForce1_m;

	// ----- solid phase -----
	double mass_s;
	double a1_s;
	double v1_s; // velocity
	double u1_s; // displacement
	
	// ----- fluid phase -----
	double volume;
	double mass_f;
	double internalForce1_f;
	double externalForce1_f;
	double nodalForce1_f;
	// w * volume = nodalForce1_f = externalForce1_f - internalForce1_f
	double w1;
	double u1_f;
};

#endif