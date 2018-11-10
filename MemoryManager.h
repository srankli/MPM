#ifndef __MEMEORYMANAGER_H__
#define __MEMEORYMANAGER_H__

/*
	The memory space returned from alloc() is raw.
	May need to be called placement new on it.
*/

class FixedSizeMemeory
{
protected:
	const size_t item_size;
	size_t cur_mem_size;
	size_t cur_pos;
	char *mem;
public:
	FixedSizeMemeory(const size_t size_of_item);
	~FixedSizeMemeory();
	inline void *get_mem() noexcept { return mem; }
	// Allocate memory with size of item_size
	void *alloc();
	// Clear memory
	void clear();
	/* Reduce the size of memory to
	get rid of excessive memory. */
	void compress();
	void reserve(size_t size);
	inline size_t capacity() noexcept { return cur_mem_size; }
	inline void *get_item_by_id(size_t id) noexcept
	{
		size_t pos_tmp = id * item_size;
		return pos_tmp < cur_pos ? mem + pos_tmp : nullptr;
	}
	inline void *get_first(void) noexcept
	{
		return cur_pos ? mem : nullptr;
	}
	inline void *get_next(void *p_item) noexcept
	{
		char *tmp = (char *)p_item + item_size;
		return tmp < mem + cur_pos ? tmp : nullptr;
	};
	inline size_t get_item_size() noexcept
	{
		return item_size;
	}
	inline size_t get_item_num(void) noexcept
	{
		return cur_pos / item_size;
	}
};


class FlexibleSizeMemory
{
protected:
	size_t cur_mem_size;
	size_t cur_pos;
	char *mem;
public:
	FlexibleSizeMemory();
	~FlexibleSizeMemory();
	// Allocate memory with size of item_size
	void *alloc(const size_t item_size);
	// Clear memory
	void clear();
	/* Reduce the size of memory to
	get rid of excessive memory. */
	void compress();
	void reserve(size_t size);
	inline size_t capacity() noexcept { return cur_mem_size; }
	inline void *get_first()
	{
		char *tmp = mem + sizeof(size_t);
		return tmp < mem + cur_pos ? tmp : nullptr;
	}
	inline void *get_next(void *p_item)
	{
		size_t si= *(size_t *)((char *)p_item - sizeof(size_t));
		char *tmp = (char *)p_item + si + sizeof(size_t);
		return tmp < mem + cur_pos ? tmp : nullptr;
	};
	inline size_t get_item_size(void *p_item) noexcept
	{
		return  *(size_t *)((char *)p_item - sizeof(size_t));
	}
};

#endif