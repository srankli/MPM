#ifndef  __ELEMENT_1D2_H__
#define  __ELEMENT_1D2_H__

#include "Mesh.h"

struct Node_1D;

/*
Local node number for each element:
	 ________________
	1                2
edge 1: node 1 ---- node 2;
*/

struct Element_1D2 : public Element
{
public:
	/*
	These two index start from 0,
	xIndex            is the position of the node1;
	xIndex + 1        is the position of the node2;
	*/
	size_t xIndex;
	/*
	Pointer to each node, these information are redundant, they are set
	for efficiency
	*/
	union
	{
		Node_1D *node[2];
		struct
		{
			Node_1D *node1;
			Node_1D *node2;
		};
	};
};

#endif