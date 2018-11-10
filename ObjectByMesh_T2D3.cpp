#include <new>

#include "ObjectByMesh_T2D3.h"

void ObjectByMesh_T2D3::addNode(Node_2D &node)
{
	Node_2D *nodeTmp;
	nodeTmp = new (static_cast<Node_2D *>(nodes_mem.alloc())) Node_2D;
	nodeTmp->index = ++curNodeIndex;
	nodeTmp->x = node.x;
	nodeTmp->y = node.y;
	++nodeNum;
}

void ObjectByMesh_T2D3::addElement(Element_T2D3 &elem)
{
	Element_T2D3 *elementTmp;
	elementTmp = new (static_cast<Element_T2D3 *>(elements_mem.alloc())) Element_T2D3;
	elementTmp->index = ++curElementIndex;
	elementTmp->node1 = elem.node1;
	elementTmp->node2 = elem.node2;
	elementTmp->node3 = elem.node3;
	++elementNum;
}
