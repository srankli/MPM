#ifndef  __ELEMENT_R2D4_H__
#define  __ELEMENT_R2D4_H__

#include "Mesh.h"

struct Node_R2D;

/*
Local node number for each element:
    4________________3
    |                |
	|                |
	|                |
	|________________|
	1                2
edge 1: node 1 ---- node 2;
edge 2: node 2 ---- node 3;
edge 3: node 3 ---- node 4;
edge 4: node 4 ---- node 1;

local coordinate:
	eta
    /\
	|
	|
	O -----> xi

*/

// 2D retangular element with 4 nodes (linear interpolation function)
struct Element_R2D4 : public Element
{
public:
	/*
	These two index start from 0,
	xIndex * mesh->nodeXNum + yIndex           is the position of the node1;
	xIndex * mesh->nodeXNum + yIndex + 1       is the position of the node2;
	(xIndex + 1) * mesh->nodeXNum + yIndex + 1 is the position of the node3;
	(xIndex + 1) * mesh->nodeXNum + yIndex     is the position of the node4.
	*/
	size_t xIndex;
	size_t yIndex;
	/*
	Pointer to each node, these information are redundant, they are set
	for efficiency
	*/
	union
	{
		Node_R2D *node[4];
		struct
		{
			Node_R2D *node1;
			Node_R2D *node2;
			Node_R2D *node3;
			Node_R2D *node4;
		};
	};
};

#endif