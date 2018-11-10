#ifndef __OBJECTBYMESH_T2D3_H__
#define __OBJECTBYMESH_T2D3_H__

#include <vector>

#include "Mesh.h"
#include "MemoryManager.h"

// Implement FEM here in the future...

class ObjectByParticle_2D_Mechanics;

class ObjectByMesh
{
protected:
	FixedSizeMemeory nodes_mem;
	size_t nodeNum;
	FixedSizeMemeory elements_mem;
	size_t elementNum;

public:
	ObjectByMesh(size_t node_size, size_t elem_size) :
		nodes_mem(node_size), nodeNum(0),
		elements_mem(elem_size), elementNum(0) {}
	~ObjectByMesh()
	{
		nodes_mem.clear();
		elements_mem.clear();
	}
};

struct Node_2D : public Node
{
	// coordinates
	union
	{
		double coords[2];
		struct { double x, y; };
	};
};

struct Node_2D_Mechanics : public Node_2D
{
	union
	{
		double velocity[2];
		struct { double v1, v2; };
	};
	Node_2D_Mechanics() : v1(0.0), v2(0.0) {}
};

struct Element_T2D3 : public Element
{
	size_t node1, node2, node3;
};

struct Element_T2D3_Mechanics : public Element_T2D3
{
	double density; // ***
	// stress, take tension as positive, default = 0.0
	union
	{
		double stress[6];
		struct
		{
			double stress11, stress22, stress33;
			double stress12, stress23, stress31;
		};
	};
	// total strain of soil skeleton, default = 0.0
	union
	{
		double strain[3];
		struct { double strain11, strain22, strain12; };
	};
	// elastic strain of soil skeleton, default = 0.0
	union
	{
		double estrain[3];
		struct { double estrain11, estrain22,  estrain12; };
	};
	// plastic strain of soil skeleton, default = 0.0
	union
	{
		double pstrain[3];
		struct { double pstrain11, pstrain22, pstrain12; };
	};
};


class ObjectByMesh_T2D3 : public ObjectByMesh
{
protected:
	size_t curNodeIndex;
	size_t curElementIndex;
public:
	ObjectByMesh_T2D3() :
		curNodeIndex(0), curElementIndex(0),
		ObjectByMesh(sizeof(Node_2D), sizeof(Element_T2D3)) {}
	~ObjectByMesh_T2D3() {}

	void addNode(Node_2D &node);
	void addElement(Element_T2D3 &elem);

	inline Node_2D *getNodeById(size_t id) noexcept
	{
		return id && id <= nodeNum ? static_cast<Node_2D *>(nodes_mem.get_mem()) + id - 1 : nullptr;
	}

	inline Element_T2D3 *getElementById(size_t id) noexcept
	{
		return id && id <= elementNum ? static_cast<Element_T2D3 *>(elements_mem.get_mem()) + id - 1 : nullptr;
	}

};


#endif