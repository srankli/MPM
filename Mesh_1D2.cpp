#include <string>
#ifdef _DEBUG
#include <iostream>
#endif

#include "Mesh_1D2.h"

#include "utilities.hpp"
#include "BoundaryCondition.h"

Mesh_1D2::~Mesh_1D2()
{
	if (nodeCoords) { delete[] nodeCoords; nodeCoords = nullptr; }
}

int Mesh_1D2::initMesh(double* nxc, size_t xnum)
{
	size_t i;
	Node_1D *node;
	Element_1D2 *element;

	nodeCoords = new double[xnum];
	memcpy(nodeCoords, nxc, sizeof(double) * xnum);
	insertionSortAsc<double>(nodeCoords, xnum);

	// initialize background mesh
	nodeNum = xnum;
	elementNum = nodeNum - 1;
	
	// initialize nodes
	nodes_mem.reserve(sizeof(Node_1D) * nodeNum);
	for (i = 0; i < nodeNum; i++)
	{
		node = static_cast<Node_1D *>(nodes_mem.alloc());
		node->index = ++curNodeIndex;
		node->xIndex = i;
	}
	nodes = static_cast<Node_1D *>(nodes_mem.get_mem());

	// initialize elements
	elements_mem.reserve(sizeof(Element_1D2) * elementNum);
	for (i = 0; i < elementNum; i++)
	{
		element = static_cast<Element_1D2 *>(elements_mem.alloc());
		element->index = ++curElementIndex;
		element->xIndex = i;
		// node1
		element->node1 = getNodeById(i + 1);
		// node2
		element->node2 = getNodeById(i + 2);
	}
	elements = static_cast<Element_1D2 *>(elements_mem.get_mem());

	return 0;
}

void Mesh_1D2::finish_init()
{
	Mesh::finish_init();
	nodes = static_cast<Node_1D *>(nodes_mem.get_mem());
	elements = static_cast<Element_1D2 *>(elements_mem.get_mem());
	
	/*// print the whole mesh
	size_t i;
	std::cout << "Node info:" << std::endl;
	for (i = 0; i < nodeNum; i++)
	{
		std::cout << "No.: " << nodes[i].index << " xIndex: "
			<< nodes[i].xIndex << std::endl;
	}
	std::cout << "Element info:" << std::endl;
	for (i = 0; i < elementNum; i++)
	{
		std::cout << "No.: " << elements[i].index << " xIndex: " 
			<< elements[i].xIndex << std::endl;
	}
	std::cout << "Velocity BC:" << std::endl;
	for (i = 0; i < velocityBCNum; i++)
	{
		std::cout << "Node_id:" << velocityBCs[i].node->index << " Value: "
		<< velocityBCs[i].v1 << std::endl;
	}
	std::cout << "Acceleration BC:" << std::endl;
	for (i = 0; i < accelerationBCNum; i++)
	{
		std::cout << "Node_id:" << accelerationBCs[i].node->index << " Value: "
			<< accelerationBCs[i].a1 << std::endl;
	}*/
}

Element_1D2 *Mesh_1D2::findInWhichElement(double x)
{
	size_t elem_x_id;

	// check if x lies in elementBuffer
	if (elementBuffer)
	{
		elem_x_id = elementBuffer->xIndex;
		if (x >= nodeCoords[elem_x_id] && x <  nodeCoords[elem_x_id + 1])
			return elementBuffer;
	}

	if (findIndex(x, &elem_x_id))
		return nullptr;
	elementBuffer = &elements[elem_x_id];

	return elementBuffer;
}

Element_1D2 *Mesh_1D2::findInWhichElement(double x, Element_1D2 *elem)
{
	size_t elem_x_id;
	long elem_id_offset;

	// first find if x lies in elem or its adjacent elements
	if (elem)
	{
		elem_x_id = elem->xIndex;
		elem_id_offset = 0;

		if (x >= nodeCoords[elem_x_id])
		{
			if (x >= nodeCoords[elem_x_id + 1])
			{
				// check if x lie outside mesh
				if (elem_x_id + 1 == elementNum)
					return nullptr;

				if (x >= nodeCoords[elem_x_id + 2])
					goto PointNotInOrNearbyElem;

				elem_id_offset++;
			}
		}
		else
		{
			// check if x lie outside mesh
			if (elem_x_id == 0)
				return nullptr;

			if (x < nodeCoords[elem_x_id - 1])
				goto PointNotInOrNearbyElem;

			elem_id_offset--;
		}

		return &elements[elem_x_id + elem_id_offset];
	}
	else if (elementBuffer) // use element buffer
	{
		elem_x_id = elementBuffer->xIndex;
		if (x >= nodeCoords[elem_x_id] && x <  nodeCoords[elem_x_id + 1])
			return elementBuffer;
		
		if (findIndex(x, &elem_x_id))
			return nullptr;
		elementBuffer = &elements[elem_x_id];
		
		return &elements[elem_x_id];
	}

PointNotInOrNearbyElem:
	if (findIndex(x, &elem_x_id))
		return nullptr;

	return &elements[elem_x_id];
}

// return index of left side of interval
// note that tail is the last index + 1
int Mesh_1D2::findIndex(double x, size_t *index)
{
	size_t head = 0, tail = nodeNum - 1;
	size_t middle;

	// check if in range
	if (x < nodeCoords[head] || x >= nodeCoords[tail])
		return -1;

	// use bisection for searching
	do
	{
		middle = (head + tail) / 2;
		if (x < nodeCoords[middle])
			tail = middle;
		else if (x > nodeCoords[middle])
			head = middle;
		else
		{
			// but this situation is rare.
			*index = middle;
			return 0;
		}
	} while (head != (tail - 1));
	*index = head;

	return 0;
}

// calculate natural coordinates and shape functions
void Mesh_1D2::calNaturalCoords(Element_1D2 *elem, double x, double *xi)
{
	double xLower, xUpper;
	double xMiddle, xHalfLength;

	xLower = nodeCoords[elem->xIndex];
	xUpper = nodeCoords[elem->xIndex + 1];

	xHalfLength = (xUpper - xLower) / 2.0;
	xMiddle = (xUpper + xLower) / 2.0;

	// xi
	*xi = (x - xMiddle) / xHalfLength;
}


void test_Mesh_1D2(void)
{
	double nodeCoords[] = { 0.0, 2.0, 1.0, 5.0, 7.0, 6.0, 3.0, 4.0, 8.0, 9.0, 10.0 };
	Mesh_1D2 mesh;
	Element_1D2 *elem;
	mesh.initMesh(nodeCoords, sizeof(nodeCoords) / sizeof(nodeCoords[0]));
	
	elem = mesh.findInWhichElement(0.5);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(0.6);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(5.0);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(6.1);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(6.5);
	std::cout << elem->xIndex << std::endl;

	elem = mesh.findInWhichElement(0.5, nullptr);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(0.6, nullptr);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(1.1, nullptr);
	std::cout << elem->xIndex << std::endl;

	elem = mesh.findInWhichElement(0.5, elem);
	std::cout << elem->xIndex << std::endl;
	elem = mesh.findInWhichElement(2.1, elem);
	std::cout << elem->xIndex << std::endl;

	double xi;
	mesh.calNaturalCoords(elem, 2.2, &xi);
	std::cout << xi << std::endl;

}
