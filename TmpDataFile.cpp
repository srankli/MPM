#include <cassert>
#ifdef _DEBUG
#include <iostream>
#endif

#include "TmpDataFile.h"

const char TmpDataFile::file_tag[] = "Numagic!";

TmpDataFile::~TmpDataFile()
{
	closeFile();
}

// mode: 0 read, 1 write
int TmpDataFile::initFile(const char *filename, unsigned char mode)
{
	assert(filename);

	size_t file_name_len = strlen(filename);
	if (file_name_len > 511)
		file_name_len = 511;
	memcpy(file_name, filename, file_name_len);
	file_name[file_name_len] = '\0';
	if (file_name_len < 8 || strcmp(file_name + file_name_len - 8, ".tmpdata"))
		memcpy(file_name + file_name_len, ".tmpdata\0", 9);

	closeFile();

	switch (mode)
	{
	case 0: // read
		pFile = fopen(file_name, "rb");
		if (!pFile) return -1;
		// skip the tag
		curFilePos = strlen(file_tag);
		fseek(pFile, (long)curFilePos, SEEK_SET);
		nextMetaDataPos = curFilePos;
		break;
	case 1:
		pFile = fopen(file_name, "wb+");
		if (!pFile)	return -1;
		// write tag
		fseek(pFile, 0, SEEK_SET);
		fwrite(file_tag, sizeof(char), strlen(file_tag), pFile);
		curFilePos = ftell(pFile);
		break;
	default:
		break;
	}

	return 0;
}


void TmpDataFile::closeFile()
{
	meta_data_write.clear();
	meta_data_read.clear();
	if (pFile)
	{
		fflush(pFile);
		fclose(pFile);
		pFile = nullptr;
	}
}


void TmpDataFile::deleteHandle(HTmpData handle)
{
	std::list<MetaDataRead>::iterator iter0;
	std::list<MetaDataWrite>::iterator iter1;
	
	switch (reinterpret_cast<MetaDataHead *>(handle)->mode)
	{
	case 0:
		iter0 = std::find(meta_data_read.begin(), meta_data_read.end(),
			*(reinterpret_cast<MetaDataRead *>(handle)));
		meta_data_read.erase(iter0);
	case 1:
		iter1 = std::find(meta_data_write.begin(), meta_data_write.end(),
			*(reinterpret_cast<MetaDataWrite *>(handle)));
		meta_data_write.erase(iter1);
	default:
		break;
	}
}

void TmpDataFile::flushFile()
{
	fflush(pFile);
}

const HTmpData TmpDataFile::addMetaData(const char *data_name,
	const void *meta_data, size_t meta_data_info_size)
{
	size_t data_type_len;
	MetaDataInFile_FixSizePart meta_data_fix_tmp;
	MetaDataWrite meta_data_tmp;

	meta_data_fix_tmp.next_meta_data = 0;
	meta_data_fix_tmp.first_data_chunk = 0;
	memset(meta_data_fix_tmp.data_type, 0, 64);
	memset(meta_data_tmp.data_type, 0, 64);
	meta_data_fix_tmp.meta_data_info_size = meta_data_info_size;
	if (data_name)
	{
		data_type_len = strlen(data_name);
		if (data_type_len > 64) data_type_len = 64;
		memcpy(meta_data_fix_tmp.data_type, data_name, data_type_len);
		memcpy(meta_data_tmp.data_type, data_name, data_type_len);
	}

	meta_data_tmp.curDataPointer = curFilePos + offsetof(MetaDataInFile_FixSizePart, first_data_chunk);

	// set pointer
	if (curMetaDataPointer)
	{
		fseek(pFile, (long)curMetaDataPointer, SEEK_SET);
		fwrite(&curFilePos, sizeof(curFilePos), 1, pFile);
	}
	curMetaDataPointer = curFilePos + offsetof(MetaDataInFile_FixSizePart, next_meta_data);

	fseek(pFile, (long)curFilePos, SEEK_SET);
	fwrite(&meta_data_fix_tmp, sizeof(meta_data_fix_tmp), 1, pFile);
	if (meta_data && meta_data_info_size)
		fwrite(meta_data, sizeof(char), meta_data_info_size, pFile);
	else
		meta_data_fix_tmp.meta_data_info_size = 0;
	
	curFilePos = ftell(pFile);

	meta_data_write.push_back(meta_data_tmp);

	return reinterpret_cast<HTmpData>(&(meta_data_write.back()));
}

void TmpDataFile::addDataChunk(const HTmpData meta_data,
	const void *data, size_t data_size)
{
	DataChunkInFile_FixSizePart data_chunk_fix_tmp;
	MetaDataWrite *pMetaData;
	size_t ori_data_size;
	size_t ori_data_size_pos;

	if (meta_data)
	{
		pMetaData = reinterpret_cast<MetaDataWrite *>(meta_data);

		if (meta_data == lastMetaDataHandle)
		{
			// update the original data chunk size
			ori_data_size_pos = pMetaData->curDataPointer
				+ offsetof(DataChunkInFile_FixSizePart, data_size);
			fseek(pFile, (long)ori_data_size_pos, SEEK_SET);
			fread(&ori_data_size, sizeof(ori_data_size), 1, pFile);
			ori_data_size += data_size;
			fseek(pFile, (long)ori_data_size_pos, SEEK_SET);
			fwrite(&ori_data_size, sizeof(ori_data_size), 1, pFile);
			
			fseek(pFile, (long)curFilePos, SEEK_SET);
			fwrite(data, sizeof(char), data_size, pFile);
			curFilePos = ftell(pFile);
		}
		else
		{
			data_chunk_fix_tmp.next_data_chunk = 0;
			data_chunk_fix_tmp.data_size = data_size;

			fseek(pFile, (long)pMetaData->curDataPointer, SEEK_SET);
			fwrite(&curFilePos, sizeof(curFilePos), 1, pFile);
			pMetaData->curDataPointer = curFilePos 
				+ offsetof(DataChunkInFile_FixSizePart, next_data_chunk);

			fseek(pFile, (long)curFilePos, SEEK_SET);
			fwrite(&data_chunk_fix_tmp, sizeof(data_chunk_fix_tmp), 1, pFile);
			fwrite(data, sizeof(char), data_size, pFile);
			curFilePos = ftell(pFile);

			lastMetaDataHandle = meta_data;
		}
	}
}

HTmpData TmpDataFile::getNextMetaData()
{
	MetaDataInFile_FixSizePart meta_data_header;
	MetaDataRead meta_data;
	
	if (!nextMetaDataPos)
		return nullptr;

	fseek(pFile, (long)nextMetaDataPos, SEEK_SET);
	fread(&meta_data_header, sizeof(meta_data_header), 1, pFile);
	nextMetaDataPos = meta_data_header.next_meta_data;

	meta_data.metaDataInfoPos = ftell(pFile);
	meta_data.metaDataInfoLen = meta_data_header.meta_data_info_size;
	meta_data.firstDataChunkPos = meta_data_header.first_data_chunk;
	meta_data.nextDataChunkPos = meta_data_header.first_data_chunk;
	memcpy(meta_data.data_type, meta_data_header.data_type, 64);
	meta_data.curDataLen = 0;

	meta_data_read.push_back(meta_data);
	return reinterpret_cast<HTmpData>(&(meta_data_read.back()));
}

const char *TmpDataFile::getDataType(HTmpData meta_data)
{
	if (meta_data) return reinterpret_cast<MetaDataHead *>(meta_data)->data_type;
	return nullptr;
}

size_t TmpDataFile::getMetaDataInfoLen(HTmpData meta_data)
{
	if (meta_data)
		return reinterpret_cast<MetaDataRead *>(meta_data)->metaDataInfoLen;
	return 0;
}

void TmpDataFile::getMetaDataInfo(HTmpData meta_data, void *meta_data_info)
{
	MetaDataRead *pMetaData;
	if (meta_data)
	{
		pMetaData = reinterpret_cast<MetaDataRead *>(meta_data);
		if (pMetaData->metaDataInfoLen)
		{
			fseek(pFile, (long)pMetaData->metaDataInfoPos, SEEK_SET);
			fread(meta_data_info, pMetaData->metaDataInfoLen, 1, pFile);
		}
	}
}

void TmpDataFile::resetDataChunk(HTmpData meta_data)
{
	MetaDataRead *pMetaData;
	
	if (meta_data)
	{
		pMetaData = reinterpret_cast<MetaDataRead *>(meta_data);
		pMetaData->nextDataChunkPos = pMetaData->firstDataChunkPos;
	}
}

size_t TmpDataFile::getNextDataChunk(HTmpData meta_data)
{
	DataChunkInFile_FixSizePart data_chunk_header;
	MetaDataRead *pMetaData;

	if (meta_data)
	{
		pMetaData = reinterpret_cast<MetaDataRead *>(meta_data);
		if (pMetaData->nextDataChunkPos)
		{
			fseek(pFile, (long)pMetaData->nextDataChunkPos, SEEK_SET);
			fread(&data_chunk_header, sizeof(data_chunk_header), 1, pFile);
			pMetaData->curDataPos = ftell(pFile);
			pMetaData->curDataLen = data_chunk_header.data_size;
			pMetaData->nextDataChunkPos = data_chunk_header.next_data_chunk;
			return pMetaData->curDataLen;
		}
	}

	return 0;
}

size_t TmpDataFile::getDataLen(HTmpData meta_data)
{
	if (meta_data)
		return reinterpret_cast<MetaDataRead *>(meta_data)->curDataLen;
	return 0;
}

void TmpDataFile::getData(HTmpData meta_data, void *data)
{
	MetaDataRead *pMetaData;
	
	if (meta_data)
	{
		pMetaData = reinterpret_cast<MetaDataRead *>(meta_data);
		if (pMetaData->curDataLen)
		{
			fseek(pFile, (long)pMetaData->curDataPos, SEEK_SET);
			fread(data, pMetaData->curDataLen, 1, pFile);
		}
	}
}


#ifdef _DEBUG

void test_ReadTmpDataFile()
{
	TmpDataFile f;
	size_t data_info_len;
	char *data_info;
	size_t data_len;
	char *data;
	HTmpData meta_data1, meta_data2, meta_data3, meta_data4;
	size_t i;

	f.initFile("ddd.tmpdata", 0);
	meta_data1 = f.getNextMetaData();
	meta_data2 = f.getNextMetaData();
	meta_data3 = f.getNextMetaData();
	meta_data4 = f.getNextMetaData();
	std::cout << meta_data3 << std::endl;
	std::cout << meta_data4 << std::endl;

	// get meta data info
	data_info_len = f.getMetaDataInfoLen(meta_data1);
	data_info = new char[data_info_len];
	f.getMetaDataInfo(meta_data1, data_info);
	for (i = 0; i < 10; i++)
		std::cout << (unsigned short)data_info[i] << " ";
	std::cout << std::endl;
	delete[] data_info;

	// get data of meta_data1
	data_len = f.getNextDataChunk(meta_data1);
	std::cout << data_len << " ";
	data_len = f.getDataLen(meta_data1);
	std::cout << data_len << std::endl;
	data = new char[data_len];
	f.getData(meta_data1, data);
	for (i = 0; i < data_len; i++)
		std::cout << data[i];
	std::cout << std::endl;
	delete[] data;

	// get data of meta_data2
	data_len = f.getNextDataChunk(meta_data2);
	std::cout << data_len << " ";
	data_len = f.getDataLen(meta_data2);
	std::cout << data_len << std::endl;
	data = new char[data_len];
	f.getData(meta_data2, data);
	for (i = 0; i < data_len; i++)
		std::cout << data[i];
	std::cout << std::endl;
	delete[] data;
	data_len = f.getNextDataChunk(meta_data2);
	std::cout << data_len << std::endl;

	f.resetDataChunk(meta_data1);
	// get data of meta_data1
	data_len = f.getNextDataChunk(meta_data1);
	std::cout << data_len << std::endl;
	data = new char[data_len];
	f.getData(meta_data1, data);
	for (i = 0; i < data_len; i++)
		std::cout << data[i];
	std::cout << std::endl;
	delete[] data;
	
	data_len = f.getNextDataChunk(meta_data1);
	std::cout << data_len << std::endl;
}


void test_WriteTmpDataFile(void)
{
	TmpDataFile f;
	char data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	char data_chunk[] = "ABCDEFG";
	HTmpData meta_data1, meta_data2, meta_data3;

	f.initFile("ddd", 1);
	meta_data1 = f.addMetaData("Test Data1", data, sizeof(data));
	meta_data2 = f.addMetaData("Test Data2", data, sizeof(data));
	f.addDataChunk(meta_data1, data_chunk, strlen(data_chunk));
	f.addDataChunk(meta_data1, data_chunk, strlen(data_chunk));
	meta_data3 = f.addMetaData("Test Data3", data, sizeof(data));
	f.addDataChunk(meta_data2, data_chunk, strlen(data_chunk));

	f.deleteHandle(meta_data1);
	f.deleteHandle(meta_data2);
	//std::cout << f.getnum() << std::endl;
}

#endif