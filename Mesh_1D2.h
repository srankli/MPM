#ifndef __MESH_1D2_H__
#define __MESH_1D2_H__

#include "Node_1D.h"
#include "Element_1D2.h"

//class Solver_1D_Mechanics_1D2_Explicit;
class Solver_1D_Mechanics_1D2_Explicit_FixedMem;
class Solver_1D_Hydromechanics_1D2_Explicit_FixedMem;

struct VelocityBCParam_1D;
struct AccelerationBCParam_1D;

class Mesh_1D2 : public Mesh
{
	//friend Solver_1D_Mechanics_1D2_Explicit;
	friend Solver_1D_Mechanics_1D2_Explicit_FixedMem;
	friend Solver_1D_Hydromechanics_1D2_Explicit_FixedMem;
protected:
	// Coordinates of nodes
	double *nodeCoords;
	Node_1D *nodes;
	Element_1D2 *elements;

public:
	Mesh_1D2() : 
		Mesh(sizeof(Node_1D), sizeof(Element_1D2), MeshType::Mesh_1D2),
		nodes(nullptr), elements(nullptr),
		nodeCoords(nullptr), elementBuffer(nullptr) {}
	~Mesh_1D2();
	// node that this Mesh object directly use data in nxc and nyc
	int initMesh(double* nxc, size_t xnum);
	void finish_init();
	// Can only be called after finish_init() is called
	inline Node_1D *getNodes(void) noexcept { return nodes; }
	inline double *getNodeCoords(void) noexcept { return nodeCoords; }
	inline Element_1D2 *getElements(void) noexcept { return elements; }

protected:
	// for accelerating inWhichElement(double *coords;)
	Element_1D2 *elementBuffer;
	// Find in which coords interval the coordinate lies
	// return -1 if coordinate lie out side the cordinates range
	int findIndex(double x, size_t *index);
public:
	// find point locate in which element;
	Element_1D2 *findInWhichElement(double x);
	// first search if point line in elem and its adjacent element
	Element_1D2 *findInWhichElement(double x, Element_1D2 *elem);

	// calculate natural coordinate
	void calNaturalCoords(Element_1D2 *elem, double x, double *xi);

	inline Node_1D *getNodeById(size_t id) noexcept
	{
		return id && id <= nodeNum ?  nodes + id - 1 : nullptr;
	}
	inline Element_1D2 *getElementById(size_t id) noexcept
	{
		return id && id <= elementNum ? elements + id - 1 : nullptr;
	}

	// ---------------------- shape function ------------------------
	inline static double N1(double *naturalCoords) noexcept
	{
		return (1.0 - naturalCoords[0]) / 2.0;
	}

	inline static double N2(double *naturalCoords) noexcept
	{
		return (1.0 + naturalCoords[0]) / 2.0;
	}

	inline static double N1(double xi) noexcept
	{
		return (1.0 - xi) / 2.0;
	}

	inline static double N2(double xi) noexcept
	{
		return (1.0 + xi) / 2.0;
	}

	// --------- one order derivative of shape function ---------
	inline static double dN1_dxi(double *naturalCoords) noexcept
	{
		return -1.0 / 2.0;
	}

	inline static double dN2_dxi(double *naturalCoords) noexcept
	{
		return 1.0 / 2.0;
	}

	inline static double dN1_dxi() noexcept
	{
		return -1.0 / 2.0;
	}

	inline static double dN2_dxi() noexcept
	{
		return 1.0 / 2.0;
	}
};

#endif