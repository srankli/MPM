#ifndef __TMPDATAFILE_H__
#define __TMPDATAFILE_H__

#include <list>
#include <cstdlib>
#include <algorithm>

// Only support file name with no longer than 512 byte

/* 
Each HtmpData point to one meta data in
TmpDataFile, can be either MetaDataRead
or MetaDataWrite
*/
typedef void* HTmpData;

class TmpDataFile
{
	// Handle of meta data
	// can be pointer of either MetaDataWrite or MetaDataRead
protected:
	struct MetaDataHead
	{
		char data_type[64];
		// 0 read, 1 write;
		const unsigned int mode;
		MetaDataHead(unsigned int md) : mode(md) {}
	};
	// for read utilities
	struct MetaDataRead
	{
		union
		{
			MetaDataHead head;
			struct
			{
				char data_type[64];
				unsigned int mode;
			};
		};
		size_t metaDataInfoPos;
		size_t metaDataInfoLen;
		size_t curDataPos;
		size_t curDataLen;
		size_t firstDataChunkPos;
		size_t nextDataChunkPos;
		bool operator == (const MetaDataRead &mdr) noexcept
		{
			return this == &mdr;
		}
		MetaDataRead() : head(0) {}
	};
	// for write utilities
	struct MetaDataWrite
	{
		union
		{
			MetaDataHead head;
			struct
			{
				char data_type[64];
				unsigned int mode;
			};
		};
		size_t curDataPointer;
		bool operator == (const MetaDataWrite &mdw) noexcept
		{
			return this == &mdw;
		}
		MetaDataWrite() : head(1) {}
	};

protected:
	// This is the tag at the head of each temporary data files
	static const char file_tag[];
	char file_name[512];
	FILE *pFile;

	// for read utilities
	std::list<MetaDataRead> meta_data_read;
	size_t nextMetaDataPos;

	// for write utilities
	std::list<MetaDataWrite> meta_data_write;
	size_t curMetaDataPointer;
	size_t curFilePos;
	HTmpData lastMetaDataHandle;

public:
	TmpDataFile() : pFile(nullptr),
		curMetaDataPointer(0), curFilePos(0), lastMetaDataHandle(0),
		nextMetaDataPos(0) {}
	~TmpDataFile();
	size_t getnum() { return meta_data_write.size(); }

public:
	struct MetaDataInFile_FixSizePart
	{
		unsigned long long next_meta_data;
		unsigned long long first_data_chunk;
		// The meaning of this tag is defined by user
		char data_type[64];
		unsigned long long meta_data_info_size;
	};
	struct DataChunkInFile_FixSizePart
	{
		unsigned long long next_data_chunk;
		unsigned long long data_size;
	};

public:
	// mode: 0 read, 1 write
	int initFile(const char *filename, unsigned char mode);
	void closeFile();
	void deleteHandle(HTmpData handle);
	const char *getFileName() { return file_name; }
	// --------------------  Read utilities --------------------
	HTmpData getNextMetaData();
	const char *getDataType(HTmpData meta_data);
	size_t getMetaDataInfoLen(HTmpData meta_data);
	void getMetaDataInfo(HTmpData meta_data, void *meta_data_info);
	void resetDataChunk(HTmpData meta_data);
	size_t getNextDataChunk(HTmpData meta_data);
	size_t getDataLen(HTmpData meta_data);
	void getData(HTmpData meta_data, void *data);

	// -------------------- Write utilities -------------------
	const HTmpData addMetaData(const char *data_name,
		const void *meta_data = nullptr, size_t meta_data_info_size = 0);
	void addDataChunk(const HTmpData meta_data, const void *data, size_t data_size);
	void flushFile();
};

#endif