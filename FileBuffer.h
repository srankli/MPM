#ifndef __FILEBUFFER_H__
#define __FILEBUFFER_H__

#include "TmpDataFile.h"

/* ------------------------- FileBuffer ------------------------- */
class FileBuffer
{
protected:
	char *buffer;
	size_t buffer_size; // = buffer_row_len * buffer_row_num
	size_t cur_pos;
	TmpDataFile *file;
	HTmpData data_handle;

public:
	FileBuffer() : buffer(nullptr), buffer_size(0),
		cur_pos(0), file(nullptr), data_handle(0) {}
	~FileBuffer() { if (buffer) delete[] buffer; }
	int initBuffer(size_t buf_size,
		const TmpDataFile *outFile, HTmpData meta_data)
	{
		if (!buf_size || !outFile || !meta_data) return -1;
		buffer = new char[buf_size];
		if (!buffer) return -2;
		buffer_size = buf_size;
		file = const_cast<TmpDataFile *>(outFile);
		data_handle = meta_data;
		return 0;
	}
	void flushBuffer()
	{
		if (cur_pos) // only flush when there is data in buffer.
			file->addDataChunk(data_handle, buffer, sizeof(char) * cur_pos);
		cur_pos = 0;
	}
#ifdef _DEBUG
	void printBuffer(void)
	{
		std::cout << "cur_pos: " << cur_pos << std::endl;
		for (size_t i = 0; i < buffer_size; i++)
			std::cout << (int)buffer[i] << " ";
		std::cout << std::endl;
	}
#endif
	// return memory of one row in the buffer so
	// that it can be filled by the caller
	size_t writeData(const void *data, size_t data_len)
	{
		size_t prev_pos, len1, len2;

		if (!data_len) return 0;

		// if data len >= buffer_size, then output data directly
		if (data_len >= buffer_size)
		{
			flushBuffer();
			//printBuffer();
			file->addDataChunk(data_handle, data, sizeof(char) * data_len);
			return data_len;
		}

		prev_pos = cur_pos;
		cur_pos += data_len;
		if (cur_pos > buffer_size)
		{
			len1 = buffer_size - prev_pos;
			len2 = cur_pos - buffer_size;
			// first fill the buffer and output to the file
			if (len1 > 0) memcpy(buffer + prev_pos, data, len1);
			//printBuffer();
			file->addDataChunk(data_handle, buffer, sizeof(char) * buffer_size);
			// add remaining data
			if (len2 > 0) memcpy(buffer, (char*)data + len1, len2);
			cur_pos = len2;
			//printBuffer();
		}
		else
		{
			memcpy(buffer + prev_pos, data, sizeof(char) * data_len);
			//printBuffer();
		}

		return data_len;
	}
};

// content len + content
class AttributeBuffer : public FileBuffer
{
public:
	size_t addAttribute(const void *data, size_t data_len)
	{
		size_t res;
		unsigned long long len_tmp;
		len_tmp = (unsigned long long)data_len;
		res = writeData(&len_tmp, sizeof(unsigned long long));
		res += writeData(data, data_len);
		return res;
	}
};

template<typename DataType>
class DataBuffer : public FileBuffer
{
protected:
	size_t buffer_row_len;
	size_t buffer_row_num;
	size_t row_byte;
public:
	int initBuffer(size_t row_len, size_t row_num,
		const TmpDataFile *outFile, HTmpData meta_data)
	{
		if (!row_len || !row_num) return -1;
		row_byte = row_len * sizeof(DataType);
		return FileBuffer::initBuffer(row_len * row_num * sizeof(DataType),
									  outFile, meta_data);
	}
	// return the head address of data buffer 
	DataType *getDataBuffer(void)
	{
		DataType *row;
		if (cur_pos >= buffer_size) flushBuffer(); // make cur_pos = 0;
		row = reinterpret_cast<DataType *>(buffer + cur_pos);
		cur_pos += row_byte;
		return row;
	}
};


#endif // !__FILEBUFFER_H__