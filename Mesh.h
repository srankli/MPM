#ifndef __MESH_H__
#define __MESH_H__

#include "MemoryManager.h"
#include "BoundaryCondition.h"

enum class MeshType : unsigned long long 
{
	InvalidMesh = 0,
	Mesh_1D2 = 1,
	Mesh_R2D4 = 2
};

class MeshTypeName
{
protected:
	static const char *type_name[];
	static const size_t type_num;
public:
	static const char *InvalidMesh; // 0
	static const char *Mesh_1D2;      // 1
	static const char *Mesh_R2D4; // 2

	static const char *getName(size_t id);
};

/* 
Note: All index start from 1. index == 0 indicates that the object is invalid.
*/

struct Node;
class Solver;

// ------------------------------ Background Mesh -----------------------------
// variables for calculation
struct NodeVar
{
	friend Node;
	friend Solver;

	size_t objectIndex;
	bool needCal;
	
	// Point to the node that this nodeVar belong
	Node *node;

protected:
	// for node variables stack of each node
	NodeVar *next_node;

	// for the whole problem
	//NodeVar *next_solver;
};

struct Node
{
	size_t index;

// --------------- NodeVar -----------------
	NodeVar *nodeVar;

// ------------- NodeVar Stack -------------
protected: 
	NodeVar *stack_top;
	size_t stack_num;
public:
	inline size_t get_num() noexcept { return stack_num; }
	inline NodeVar *get_top(void) noexcept { return stack_top; }
	inline void push(NodeVar *item) noexcept
	{
		item->next_node = stack_top;
		stack_top = item;
		++stack_num;
	}
	inline void reset(void) noexcept
	{
		stack_top = nullptr;
		stack_num = 0;
	}
	inline void start(NodeVar **cur) noexcept { *cur = stack_top; }
	inline void next(NodeVar **cur) noexcept { *cur = (*cur)->next_node; }
};

struct Element
{
	size_t index;
};


class Mesh
{
	friend Solver;
protected:
	// Mesh type
	MeshType mtype;
	// Nodes
	size_t curNodeIndex;
	size_t nodeNum;
	FixedSizeMemeory nodes_mem;
	// Elements
	size_t curElementIndex;
	size_t elementNum;
	FixedSizeMemeory elements_mem;
	// Boundary Conditions
	// velocity boundary condtions
	size_t velocityBCNum;
	VelocityBCManager velocityBCs_mem;
	// acceleration boundary conditions
	size_t accelerationBCNum;
	AccelerationBCManager accelerationBCs_mem;

public:
	Mesh(size_t node_size, size_t elem_size, MeshType mtp) :
		mtype(mtp),
		curNodeIndex(0), nodeNum(0), nodes_mem(node_size),
		curElementIndex(0), elementNum(0), elements_mem(elem_size),
		velocityBCNum(0), accelerationBCNum(0) {}
	~Mesh()
	{
		nodes_mem.clear();
		elements_mem.clear();
		velocityBCs_mem.clear();
		accelerationBCs_mem.clear();
	}
	inline MeshType getMeshType(void) noexcept { return mtype; }
	inline size_t getNodeNum(void) noexcept { return nodeNum; }
	inline size_t getElementNum(void) noexcept { return elementNum; }

	// Maybe needed to be redefined
	virtual bool validateNodeId(size_t id)
	{
		return id && id <= nodeNum ? true : false;
	}
	// Add Boundary conditions
	int addVelocityBC(VelocityBCParam *vbc_param)
	{
		int res = -1; 
		if (validateNodeId(vbc_param->index))
		{
			res = velocityBCs_mem.add_bc(vbc_param);
			if (!res) ++velocityBCNum;
		}
		return res;
	}
	int addAccelerationBC(AccelerationBCParam *abc_param)
	{
		int res = -1;
		if (validateNodeId(abc_param->index))
		{
			res = accelerationBCs_mem.add_bc(abc_param);
			if (!res) ++accelerationBCNum;
		}
		return res;
	}

	void finish_init()
	{
		nodes_mem.compress();
		elements_mem.compress();
		velocityBCs_mem.compress();
		accelerationBCs_mem.compress();
	}
};

#endif