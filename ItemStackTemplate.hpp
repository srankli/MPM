#ifndef __ITEMSTACKTEMPLATE_HPP__
#define __ITEMSTACKTTMPLATE_HPP__

#include "MemoryPoolTemplate.hpp"

/*
1. This stack use memory pool to enhance speed of memeory operation
2. Item should have default constructor function
*/
template <typename Item>
class ItemStackTemplate
{
public:
	struct StackNode
	{
		Item item;
		// List chain starts from the top of the stack
		// top -> StackNode_n -> StackNode_n-1 -> ...
		StackNode *prev;
		StackNode(StackNode *prevNode = nullptr)
			: prev(prevNode) {}
	};

protected:
	MemoryPoolTemplate<StackNode> memPool;

	StackNode *top; // stack top
	unsigned int num;
	StackNode *cur, *cur_0;

public:
	ItemStackTemplate(unsigned int defaultSlotSize = 0) :
		top(nullptr), num(0),
		cur(nullptr), cur_0(nullptr),
		memPool(defaultSlotSize) {}
	~ItemStackTemplate() { reset(); }
	void setSlotSize(unsigned int size)
	{
		memPool.setSlotSize(size);
	}
	// Clear the whole stack;
	void reset(void) noexcept;
	inline unsigned int getNum(void) { return num; }
	// Stack operations
	inline Item *push(void);
	inline Item *push(Item &item);
	// return -1 if there is not object to be poped
	int pop(void) noexcept; // Dump item directly
	int pop(Item &item) noexcept;
	// List operations
	inline Item *add(void); // = push()
	inline Item *getTop(void) { return top; }

	// Traversal from top of the stack
	// During Traversal operation, push, pop, add are not allowed.
	inline Item *start(void)
	{
		cur_0 = nullptr;
		cur = top;
		return reinterpret_cast<Item *>(cur);
	}
	inline Item *next(void)
	{
		cur_0 = cur;
		cur = cur ? cur->prev : nullptr;
		return reinterpret_cast<Item *>(cur);
	}
	// return -1 if fail, 0 if success
	int delCurItem(void) noexcept;

#ifdef _DEBUG
	void printStack(void);
#endif
};

template <typename Item>
void ItemStackTemplate<Item>::reset(void) noexcept
{
	/*
	StackNode *tmp;
	tmp = top;
	while (tmp)
	{
		top = top->prev;
		memPool->free(tmp);
		tmp = top;
	}*/
	memPool.clear();
	top = nullptr;
	num = 0;
	cur = nullptr;
	cur_0 = nullptr;
}

/*
template <typename Item>
void ItemStackTemplate<Item>::reset_notfree(void) noexcept
{
	top = nullptr;
	num = 0;
	cur = nullptr;
	cur_0 = nullptr;
}*/

template <typename Item>
Item *ItemStackTemplate<Item>::push(void)
{
	StackNode *tmp;
	++num;
	tmp = memPool.alloc();
	tmp->prev = top;
	top = tmp;
	return &(tmp->item);
}

template <typename Item>
Item *ItemStackTemplate<Item>::push(Item &item)
{
	StackNode *tmp;
	++num;
	tmp = memPool.alloc();
	tmp->item = item;
	tmp->prev = top;
	top = tmp;
	return &(tmp->item);
}

template <typename Item>
int ItemStackTemplate<Item>::pop() noexcept
{
	StackNode *tmp;
	if (top)
	{
		--num;
		tmp = top;
		top = top->prev;
		memPool.free(tmp);
		return 0;
	}
	return -1;
}

template <typename Item>
int ItemStackTemplate<Item>::pop(Item &item) noexcept
{
	StackNode *tmp;
	if (top)
	{
		--num;
		/*
		if (cur == top)
			cur = top->prev;
		if (cur_0 == top)
		{
			cur_0 = top->prev;
			cur = cur ? cur->prev : nullptr;
		}*/
		item = top->item;
		tmp = top;
		top = top->prev;
		memPool.free(tmp);
		return 0;
	}
	return -1;
}

template <typename Item>
Item *ItemStackTemplate<Item>::add(void) // = push()
{
	StackNode *tmp;
	++num;
	tmp = memPool.alloc();
	tmp->prev = top;
	top = tmp;
	return &(tmp->item);
}

// The cost of del is relatively high
template <typename Item>
int ItemStackTemplate<Item>::delCurItem(void) noexcept
{
	if (cur)
	{
		// cur does not point to the top item
		if (cur_0)
		{
			cur_0->prev = cur->prev;
			memPool.free(cur);
			cur = cur_0->prev;
		}
		else // cur points to the first item
		{
			top = cur->prev;
			memPool.free(cur);
			cur = top;
		}
		return 0;
	}
	return -1;
}

#ifdef _DEBUG
template <typename Item>
void ItemStackTemplate<Item>::printStack(void)
{
	StackNode *tmp;
	for (tmp = top; tmp; tmp = tmp->prev)
	{
		std::cout << tmp->item << " ";
	}
	std::cout << std::endl;
}
#endif

#endif