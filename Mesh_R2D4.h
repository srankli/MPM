#ifndef __MESH_R2D_H__
#define __MESH_R2D_H__

#include "Node_R2D.h"
#include "Element_R2D4.h"

class Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem;

struct VelocityBCParam_2D;
struct AccelerationBCParam_2D;

extern int main(void);

class Mesh_R2D4 : public Mesh
{
	friend Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem;
protected:
	// Coordinates of nodes
	double *nodeXCoords;
	double *nodeYCoords;
	size_t nodeXNum;
	size_t nodeYNum;
	size_t elementXNum;
	size_t elementYNum;
	Node_R2D *nodes;
	Element_R2D4 *elements;

public:
	Mesh_R2D4() :
		Mesh(sizeof(Node_R2D), sizeof(Element_R2D4), MeshType::Mesh_R2D4),
		nodeXCoords(nullptr), nodeYCoords(nullptr),
		nodeXNum(0), nodeYNum(0),
		elementXNum(0), elementYNum(0),
		nodes(nullptr), elements(nullptr),
		elementBuffer(nullptr) {}
	~Mesh_R2D4();
	// node that this Mesh object directly use data in nxc and nyc
	int initMesh(double *nxc, size_t xnum, double *nyc, size_t ynum);
	void finish_init();
	// Can only be called after finish_init() is called
	inline Node_R2D *getNodes(void) noexcept { return nodes; }
	inline Element_R2D4 *getElements(void) noexcept { return elements; }
	inline size_t getXCoordNum(void) noexcept { return nodeXNum; }
	inline size_t getYCoordNum(void) noexcept { return nodeYNum; }
	inline double *getXNodeCoords(void) noexcept { return nodeXCoords; }
	inline double *getYNodeCoords(void) noexcept { return nodeYCoords; }

protected:
	// for accelerating inWhichElement(double *coords;)
	Element_R2D4 *elementBuffer;
	// Find in which coords interval the coordinate lies
	// return -1 if coordinate lie out side the cordinates range
	int findXIndex(double x, size_t *index);
	int findYIndex(double y, size_t *index);

public:
	// find point locate in which element;
	Element_R2D4 *findInWhichElement(double x, double y);
	// first search if point line in elem and its adjacent element
	Element_R2D4 *findInWhichElement(double x, double y, Element_R2D4 *elem);

	// calculate natural coordinate
	void calNaturalCoords(Element_R2D4 *elem, double x, double y, double *xi, double *eta);

	inline Node_R2D *getNodeById(size_t id) noexcept
	{
		return id && id <= nodeNum ? nodes + id - 1 : nullptr;
	}

	inline Element_R2D4 *getElementById(size_t id) noexcept
	{
		return id && id <= elementNum ? elements + id - 1 : nullptr;
	}

	inline Node_R2D *getNodeByXYIndex(size_t x_id, size_t y_id) noexcept
	{
		return nodes + y_id * nodeXNum + x_id;
	}

	inline Element_R2D4 *getElementByXYId(size_t x_id, size_t y_id) noexcept
	{
		return elements + y_id * elementXNum + x_id;
	}

	// ----------------------- shape function -----------------------
	inline static double N1(double *naturalCoords) noexcept
	{
		return (1.0 - naturalCoords[0]) * (1.0 - naturalCoords[1]) / 4.0;
	}

	inline static double N2(double *naturalCoords) noexcept
	{
		return (1.0 + naturalCoords[0]) * (1.0 - naturalCoords[1]) / 4.0;
	}

	inline static double N3(double *naturalCoords) noexcept
	{
		return (1.0 + naturalCoords[0]) * (1.0 + naturalCoords[1]) / 4.0;
	}

	inline static double N4(double *naturalCoords) noexcept
	{
		return (1.0 - naturalCoords[0]) * (1.0 + naturalCoords[1]) / 4.0;
	}

	inline static double N1(double xi, double eta) noexcept
	{
		return (1.0 - xi) * (1.0 - eta) / 4.0;
	}

	inline static double N2(double xi, double eta) noexcept
	{
		return (1.0 + xi) * (1.0 - eta) / 4.0;
	}

	inline static double N3(double xi, double eta) noexcept
	{
		return (1.0 + xi) * (1.0 + eta) / 4.0;
	}

	inline static double N4(double xi, double eta) noexcept
	{
		return (1.0 - xi) * (1.0 + eta) / 4.0;
	}

	// --------- one order derivative of shape function ---------
	inline static double dN1_dxi(double *naturalCoords)  noexcept { return -(1.0 - naturalCoords[1]) / 4.0; }
	inline static double dN1_deta(double *naturalCoords) noexcept {	return -(1.0 - naturalCoords[0]) / 4.0; }
	inline static double dN2_dxi(double *naturalCoords)  noexcept { return  (1.0 - naturalCoords[1]) / 4.0; }
	inline static double dN2_deta(double *naturalCoords) noexcept {	return -(1.0 + naturalCoords[0]) / 4.0; }
	inline static double dN3_dxi(double *naturalCoords)  noexcept { return  (1.0 + naturalCoords[1]) / 4.0; }
	inline static double dN3_deta(double *naturalCoords) noexcept { return  (1.0 + naturalCoords[0]) / 4.0; }
	inline static double dN4_dxi(double *naturalCoords)  noexcept { return -(1.0 + naturalCoords[1]) / 4.0; }
	inline static double dN4_deta(double *naturalCoords) noexcept { return  (1.0 - naturalCoords[0]) / 4.0; }

	inline static double dN1_dxi(double xi, double eta)  noexcept { return -(1.0 - eta) / 4.0; }
	inline static double dN1_deta(double xi, double eta) noexcept { return -(1.0 - xi)  / 4.0; }
	inline static double dN2_dxi(double xi, double eta)  noexcept { return  (1.0 - eta) / 4.0; }
	inline static double dN2_deta(double xi, double eta) noexcept { return -(1.0 + xi)  / 4.0; }
	inline static double dN3_dxi(double xi, double eta)  noexcept { return  (1.0 + eta) / 4.0; }
	inline static double dN3_deta(double xi, double eta) noexcept { return  (1.0 + xi)  / 4.0; }
	inline static double dN4_dxi(double xi, double eta)  noexcept { return -(1.0 + eta) / 4.0; }
	inline static double dN4_deta(double xi, double eta) noexcept { return  (1.0 - xi) / 4.0; }

#ifdef _DEBUG
	friend void test_mesh(void);
#endif
};

#endif