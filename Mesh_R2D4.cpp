#include <string>
#ifdef _DEBUG
#include <iostream>
#endif

#include "Mesh_R2D4.h"

#include "utilities.hpp"
#include "BoundaryCondition.h"

Mesh_R2D4::~Mesh_R2D4()
{
	if (nodeXCoords) { delete[] nodeXCoords; nodeXCoords = nullptr; }
	if (nodeYCoords) { delete[] nodeYCoords; nodeXCoords = nullptr; }
}

int Mesh_R2D4::initMesh(double *nxc, size_t xnum, double *nyc, size_t ynum)
{
	size_t i, j;
	size_t node_id_tmp;
	Node_R2D *nodeTmp;
	Element_R2D4 *elementTmp;

	nodeXNum = xnum;
	nodeXCoords = new double[xnum];
	memcpy(nodeXCoords, nxc, sizeof(double) * xnum);
	insertionSortAsc<double>(nodeXCoords, nodeXNum);

	nodeYNum = ynum;
	nodeYCoords = new double[ynum];
	memcpy(nodeYCoords, nyc, sizeof(double) * ynum);
	insertionSortAsc<double>(nodeYCoords, nodeYNum);

	// Initialize background mesh
	nodeNum = nodeXNum * nodeYNum;
	nodes_mem.reserve(sizeof(Node_R2D) * nodeNum);
	for (j = 0; j < nodeYNum; j++)
	{
		for (i = 0; i < nodeXNum; i++)
		{
			nodeTmp = static_cast<Node_R2D *>(nodes_mem.alloc());
			nodeTmp->index = ++curNodeIndex;
			nodeTmp->xIndex = i;
			nodeTmp->yIndex = j;
		}
	}
	nodes = static_cast<Node_R2D *>(nodes_mem.get_mem());

	elementXNum = nodeXNum - 1;
	elementYNum = nodeYNum - 1;
	elementNum = elementXNum * elementYNum;
	for (j = 0; j < elementYNum; j++)
	{
		for (i = 0; i < elementXNum; i++)
		{
			elementTmp = static_cast<Element_R2D4 *>(elements_mem.alloc());
			elementTmp->index = ++curElementIndex;
			elementTmp->xIndex = i;
			elementTmp->yIndex = j;
			// node 1
			node_id_tmp = j * nodeXNum + i;
			elementTmp->node1 = nodes + node_id_tmp;
			// node 2
			node_id_tmp++;
			elementTmp->node2 = nodes + node_id_tmp;
			// node3
			node_id_tmp += nodeXNum;
			elementTmp->node3 = nodes + node_id_tmp;
			// node4
			node_id_tmp--;
			elementTmp->node4 = nodes + node_id_tmp;
		}
	}
	elements = static_cast<Element_R2D4 *>(elements_mem.get_mem());

	return 0;
}

void Mesh_R2D4::finish_init()
{
	Mesh::finish_init();
	nodes = static_cast<Node_R2D *>(nodes_mem.get_mem());
	elements = static_cast<Element_R2D4 *>(elements_mem.get_mem());
}


Element_R2D4 *Mesh_R2D4::findInWhichElement(double x, double y)
{
	size_t elem_x_id, elem_y_id;

	// check if coords lies in elementBuffer
	if (elementBuffer)
	{
		elem_x_id = elementBuffer->xIndex;
		elem_y_id = elementBuffer->yIndex;
		
		if (x >= nodeXCoords[elem_x_id] &&
			x <  nodeXCoords[elem_x_id + 1] &&
			y >= nodeYCoords[elem_y_id] &&
			y <  nodeYCoords[elem_y_id + 1])
			return elementBuffer;
	}

	if (findXIndex(x, &elem_x_id) == -1
		|| findYIndex(y, &elem_y_id) == -1)
		return nullptr;
	elementBuffer = getElementByXYId(elem_x_id, elem_y_id);
	
	return elementBuffer;
}

Element_R2D4 *Mesh_R2D4::findInWhichElement(double x, double y, Element_R2D4 *elem)
{
	size_t elem_x_id, elem_y_id;
	long int elem_id_offset;

	// first check if coords is in elem and its adjacent elements
	if (elem)
	{
		elem_x_id = elem->xIndex;
		elem_y_id = elem->yIndex;
		elem_id_offset = 0;

		// x direction
		if (x >= nodeXCoords[elem_x_id])
		{
			if (x >= nodeXCoords[elem_x_id + 1])
			{
				// check if coords lie outside mesh
				if (elem_x_id + 1 == elementXNum)
					return nullptr;

				if (x >= nodeXCoords[elem_x_id + 2])
					goto PointNotInOrNearElem;

				elem_id_offset++;
			}
		}
		else
		{
			// check if coords lie outside mesh
			if (elem_x_id == 0)
				return nullptr;

			if (x < nodeXCoords[elem_x_id - 1])
				goto PointNotInOrNearElem;

			elem_id_offset--;
		}

		// y direction
		if (y >= nodeYCoords[elem_y_id])
		{
			if (y >= nodeYCoords[elem_y_id + 1])
			{
				// check if coords lie outside mesh
				if (elem_y_id + 1 == elementYNum)
					return nullptr;

				if (y >= nodeYCoords[elem_y_id + 2])
					goto PointNotInOrNearElem;

				elem_id_offset += (long int)elementXNum;
			}
		}
		else
		{
			// check if coords lie outside mesh
			if (elem_y_id == 0)
				return nullptr;

			if (y < nodeYCoords[elem_y_id - 1])
				goto PointNotInOrNearElem;

			elem_id_offset -= (long int)elementXNum;
		}
		
		return &elements[elem_y_id * elementXNum + elem_x_id + elem_id_offset];
	}

PointNotInOrNearElem:
	return findInWhichElement(x, y);
}

// return index of left side of interval
// note that tail is the last index + 1
int Mesh_R2D4::findXIndex(double x, size_t *index)
{
	size_t head = 0, tail = nodeXNum - 1;
	size_t middle;

	// check if in range
	if (x < nodeXCoords[head] || x >= nodeXCoords[tail])
		return -1;

	// use bisection for searching
	do
	{
		middle = (head + tail) / 2;
		if (x < nodeXCoords[middle])
			tail = middle;
		else if (x > nodeXCoords[middle])
			head = middle;
		else
		{
			// but this situation is rare.
			*index = middle;
			return 0;
		}
	}
	while (head != (tail - 1));
	*index = head;

	return 0;
}

// return index of left side of interval
// note that tail is the last index + 1
int Mesh_R2D4::findYIndex(double y, size_t *index)
{
	size_t head = 0, tail = nodeYNum - 1;
	size_t middle;

	// check if in range
	if (y < nodeYCoords[head] || y >= nodeYCoords[tail])
		return -1;

	// use bisection for searching
	do
	{
		middle = (head + tail) / 2;
		if (y < nodeYCoords[middle])
			tail = middle;
		else if (y > nodeYCoords[middle])
			head = middle;
		else
		{
			// but this situation is rare.
			*index = middle;
			return 0;
		}
	}
	while (head != (tail - 1));
	*index = head;

	return 0;
}


// Calculate natural coordinate
void Mesh_R2D4::calNaturalCoords(Element_R2D4 *elem,
	double x, double y, double *xi, double *eta)
{
	double xLower, xUpper, yLower, yUpper;
	double xMiddle, xHalfLength, yMiddle, yHalfLength;

	xLower = nodeXCoords[elem->xIndex];
	xUpper = nodeXCoords[elem->xIndex + 1];
	yLower = nodeYCoords[elem->yIndex];
	yUpper = nodeYCoords[elem->yIndex + 1];

	xHalfLength = (xUpper - xLower) / 2.0;
	yHalfLength = (yUpper - yLower) / 2.0;
	xMiddle = (xUpper + xLower) / 2.0;
	yMiddle = (yUpper + yLower) / 2.0;

	*xi = (x - xMiddle) / xHalfLength;
	*eta = (y - yMiddle) / yHalfLength;
}

#ifdef _DEBUG
void test_mesh(void)
{
	Mesh_R2D4 mesh;
	Element_R2D4 *elems, *elem;
	Node_R2D *nodes;
	size_t i;

	double x_coords[] = { 0.0, 1.0, 2.0, 3.0 };
	double y_coords[] = { 0.0, 1.0, 2.0, 3.0 };
	mesh.initMesh(x_coords, sizeof(x_coords) / sizeof(x_coords[0]),
		y_coords, sizeof(y_coords) / sizeof(y_coords[0]));
	elems = mesh.elements;
	nodes = mesh.nodes;

	for (i = 0; i < mesh.nodeNum; i++)
		std::cout << "n_id: " << mesh.nodes[i].index
		<< " x: " << mesh.nodeXCoords[mesh.nodes[i].xIndex]
		<< " y: " << mesh.nodeYCoords[mesh.nodes[i].yIndex] << std::endl;

	std::cout << std::endl;

	for (i = 0; i < mesh.elementNum; i++)
		std::cout << "e_id: " << mesh.elements[i].index
		<< " n1_id: " << mesh.elements[i].node1->index
		<< " n2_id: " << mesh.elements[i].node2->index
		<< " n3_id: " << mesh.elements[i].node3->index
		<< " n4_id: " << mesh.elements[i].node4->index
		<< std::endl;

	elem = mesh.findInWhichElement(0.5, 2.5, nullptr);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(1.5, 1.5, elems + 4);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(0.5, 1.5, elems + 4);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(0.5, 2.5, elems + 4);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(0.5, 3.5, elems + 4);
	if (!elem) std::cout << "null" << std::endl;
	else std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(1.5, 0.5, elems + 6);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(0.5, 2.5, elems + 3);
	std::cout << elem->index << std::endl;

	elem = mesh.findInWhichElement(-0.5, 3.5, elems + 1);
	if (!elem) std::cout << "null" << std::endl;
	else std::cout << elem->index << std::endl;

	double xi, eta;
	mesh.calNaturalCoords(elems + 4, 1.0, 1.1, &xi, &eta);
	std::cout << "xi: " << xi << " eta: " << eta << std::endl;
}
#endif