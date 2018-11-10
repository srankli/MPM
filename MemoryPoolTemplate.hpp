#ifndef __MEMEORYPOOLTEMPLATE_HPP__
#define __MEMEORYPOOLTEMPLATE_HPP__

#include <cassert>

/*
Assumption:
	Item should have default constructor function.
*/
template <typename Item>
class MemoryPoolTemplate
{
protected:
	/*
	    Each PoolItem in the pool add a "next" pointer after the
	original Item type.
	    "next" points to:
			1. the freeItemListHead when the PoolItem is in used
			2, next empty PoolItem if it is not in used
	*/
	struct PoolItem
	{
		Item item;
		PoolItem *next;
	};
	/*
	Each memory pool consists of several
	memory slots
	*/
	struct MemorySlotNode
	{
	public:
		PoolItem *pool;
		size_t size;
		MemorySlotNode *next;
	public:
		MemorySlotNode(size_t slotSize,
			MemorySlotNode *nt = nullptr) :
			size(slotSize), next(nt)
		{
			assert(size > 0);
			pool = new PoolItem[size];
		}
		~MemorySlotNode() { if (pool) delete[] pool; }
	};

protected:
	size_t newSlotSize;
	MemorySlotNode *slotListHead;

protected: 
	// The memory slot can be allocated
	MemorySlotNode *curSlot;
	PoolItem *curPool;
	size_t curPoolSize, curPoolIndex;
	/* 
	List of empty PoolItem in each slots
	Note that the address of freeItemListHeadItem itself serve as a tag: 
	  For each non-empty PoolItem, its "next" component points to
	  freeItemListHead to prove that it really lies in one of the slot of this pool.
	*/
	union
	{
		PoolItem freeItemListHeadItem;
		struct
		{
			Item freeItem;
			//List of empty Item in each slots
			PoolItem *freeItemListHead;
		};
	};

public:
	MemoryPoolTemplate(size_t size = 0);
	~MemoryPoolTemplate();
	void setSlotSize(size_t size);
	Item *alloc(void);
	void free(Item *item);
	void clear(void); // free all item but do not free memory
	void clearMem(void);
	// delete slot that are not necessary
	// void compressMem();

#ifdef _DEBUG
	void printSlotInfo(void);
#endif
};

template <typename Item>
MemoryPoolTemplate<Item>::MemoryPoolTemplate(size_t size) :
	newSlotSize(size), slotListHead(nullptr),
	freeItemListHead(nullptr),
	curSlot(nullptr), curPool(nullptr),
	curPoolSize(0), curPoolIndex(0) {}

template <typename Item>
MemoryPoolTemplate<Item>::~MemoryPoolTemplate() { clearMem(); }

template <typename Item>
void MemoryPoolTemplate<Item>::setSlotSize(size_t size)
{
	assert(size);
	newSlotSize = size;
}

template <typename Item>
void MemoryPoolTemplate<Item>::clearMem(void)
{
	MemorySlotNode *curMemSlot;
	
	while (slotListHead)
	{
		curMemSlot = slotListHead;
		slotListHead = slotListHead->next;
		delete curMemSlot;
	}

	freeItemListHead = nullptr;
	curSlot = nullptr;
	curPool = nullptr;
	curPoolSize = 0;
	curPoolIndex = 0;
}

template <typename Item>
void MemoryPoolTemplate<Item>::clear()
{
	freeItemListHead = nullptr;
	curSlot = slotListHead;
	if (curSlot)
	{
		curPool = curSlot->pool;
		curPoolSize = curSlot->size;
	}
	else
	{	// if not initialized
		curPool = nullptr;
		curPoolSize = 0;
	}
	curPoolIndex = 0;
}

template <typename Item>
Item *MemoryPoolTemplate<Item>::alloc(void)
{
	PoolItem *pPoolItem;
	Item *pItem;

	// The new SlotSize need to be > 0 in order to allocate memory
	assert(newSlotSize);

	if (freeItemListHead)
	{
		pPoolItem = freeItemListHead;
		freeItemListHead = freeItemListHead->next;
		pItem = &(pPoolItem->item);
		pPoolItem->next = &freeItemListHeadItem;
	}
	else if (curPoolIndex < curPoolSize)
	{
		pItem = &(curPool[curPoolIndex].item);
		curPool[curPoolIndex].next = &freeItemListHeadItem;
		++curPoolIndex;
	}
	else if (!curSlot) // if no memory has been allocated
	{
		assert(newSlotSize > 0);
		slotListHead = new MemorySlotNode(newSlotSize, nullptr);
		curSlot = slotListHead;
		curPool = curSlot->pool;
		curPoolSize = curSlot->size;
		curPoolIndex = 1;
		pItem = &(curPool[0].item);
		curPool[0].next = &freeItemListHeadItem;
	}
	else
	{
		assert(newSlotSize > 0);
		if (!(curSlot->next)) // all memory has been used up
			curSlot->next = new MemorySlotNode(newSlotSize, nullptr);
		curSlot = curSlot->next;
		curPool = curSlot->pool;
		curPoolSize = curSlot->size;
		curPoolIndex = 1;
		pItem = &(curPool[0].item);
		curPool[0].next = &freeItemListHeadItem;
	}
	return pItem;
}

template <typename Item>
void MemoryPoolTemplate<Item>::free(Item *item)
{
	if (item)
	{
		PoolItem *poolItemTmp = reinterpret_cast<PoolItem *>(item);
		// check if the item lie in this memory pool
		if (poolItemTmp->next == &freeItemListHeadItem)
		{
			poolItemTmp->next = freeItemListHead;
			freeItemListHead = poolItemTmp;
		}
	}
}

#ifdef _DEBUG
template <typename Item>
void MemoryPoolTemplate<Item>::printSlotInfo(void)
{
	MemorySlotNode *curNode;
	unsigned int c = 0;
	for (curNode = slotListHead; curNode; curNode = curNode->next)
	{
		++c;
		std::cout << "slotNode " << c << " size: " << curNode->size << std::endl;
	}
}
#endif

#endif