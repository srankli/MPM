#include <cstdlib>
#include <cassert>
#ifdef _DEBUG
#include <iostream>
#endif

#include "MemoryManager.h"

/* ------------------ FixedSizeMemeory ----------------- */
FixedSizeMemeory::FixedSizeMemeory(const size_t size_of_item) :
	item_size(size_of_item),
	cur_mem_size(0), cur_pos(0), mem(nullptr)
{
	assert(size_of_item);
}

FixedSizeMemeory::~FixedSizeMemeory()
{
	if (mem) free(mem);
}

// Allocate memory with size of item_size
void *FixedSizeMemeory::alloc()
{
	size_t tmp_pos;
	tmp_pos = cur_pos;
	cur_pos += item_size;
	if (cur_pos >= cur_mem_size)
	{
		// The capacity expands by 2 fold
		cur_mem_size = cur_pos + cur_pos;
		mem = (char *)realloc(mem, cur_mem_size);
		if (!mem)
		{
			cur_mem_size = 0;
			cur_pos = 0;
			return nullptr;
		}
	}
	return mem + tmp_pos;
}

// Clear memory
void FixedSizeMemeory::clear()
{
	if (mem)
	{
		free(mem);
		mem = nullptr;
		cur_mem_size = 0;
		cur_pos = 0;
	}
}

/* Reduce the size of memory to
get rid of excessive memory. */
void FixedSizeMemeory::compress()
{
	if (cur_pos != cur_mem_size)
	{
		cur_mem_size = cur_pos;
		mem = (char *)realloc(mem, cur_mem_size);
		if (!mem)
		{
			cur_mem_size = 0;
			cur_pos = 0;
		}
	}
}

void FixedSizeMemeory::reserve(size_t size)
{
	assert(size);
	if (size > cur_mem_size)
	{
		mem = (char *)realloc(mem, size);
		if (mem)
			cur_mem_size = size;
		else
		{
			cur_mem_size = 0;
			cur_pos = 0;
		}
	}
}

/* ---------------- FlexiableSizeMemory --------------- */
FlexibleSizeMemory::FlexibleSizeMemory() :
	cur_mem_size(0), cur_pos(0), mem(nullptr) {}

FlexibleSizeMemory::~FlexibleSizeMemory()
{
	if (mem) free(mem);
}

void *FlexibleSizeMemory::alloc(const size_t item_size)
{
	size_t tmp_pos;

	assert(item_size);
	tmp_pos = cur_pos;
	cur_pos += sizeof(size_t) + item_size;
	if (cur_pos >= cur_mem_size)
	{
		cur_mem_size = cur_pos + cur_pos;
		mem = (char *)realloc(mem, cur_mem_size);
		if (!mem)
		{
			cur_mem_size = 0;
			cur_pos = 0;
			return nullptr;
		}
	}
	*(size_t *)(mem + tmp_pos) = item_size;
	return mem + tmp_pos + sizeof(size_t);
}

void FlexibleSizeMemory::clear()
{
	if (mem)
	{
		free(mem);
		mem = nullptr;
		cur_mem_size = 0;
		cur_pos = 0;
	}
}

void FlexibleSizeMemory::compress()
{
	if (cur_pos != cur_mem_size)
	{
		cur_mem_size = cur_pos;
		mem = (char *)realloc(mem, cur_mem_size);
		if (!mem)
		{
			cur_mem_size = 0;
			cur_pos = 0;
		}
	}
}

void FlexibleSizeMemory::reserve(size_t size)
{
	assert(size);
	if (size > cur_mem_size)
	{
		mem = (char *)realloc(mem, size);
		if (mem)
			cur_mem_size = size;
		else
		{
			cur_mem_size = 0;
			cur_pos = 0;
		}
	}
}
