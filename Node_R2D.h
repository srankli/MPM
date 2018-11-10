#ifndef __NODE_R2D_H__
#define __NODE_R2D_H__

#include "Mesh.h"

struct NodeVar;

// node for 2D retangular element
struct Node_R2D : public Node
{
public:
	size_t xIndex;
	size_t yIndex;
};

// variables for calculation
struct NodeVar_2D_Mechanics : NodeVar
{
public:
	// Physical properties
	double mass;
	union // linear momentum
	{
		double momentum[2];
		struct
		{
			double momentum1;
			double momentum2;
		};
	};
	union // increment of linear momentum
	{
		double dMomentum[2];
		struct
		{
			double dMomentum1;
			double dMomentum2;
		};
	};
	union // velocity
	{
		double velocity[2];
		struct
		{
			double v1;
			double v2;
		};
	};
	union // displacement
	{
		double displacement[2];
		struct
		{
			double u1;
			double u2;
		};
	};
	union // internal force
	{
		double internalForce[2];
		struct
		{
			double internalForce1;
			double internalForce2;
		};
	};
	union // external force
	{
		double externalForce[2];
		struct
		{
			double externalForce1;
			double externalForce2;
		};
	};
	union // contact force
	{
		double contactForce[2];
		struct
		{
			double contactForce1;
			double contactForce2;
		};
	};
	union //nodal force
	{
		double nodalForce[2];
		struct
		{
			double nodalForce1;
			double nodalForce2;
		};
	};
};

#endif