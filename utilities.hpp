#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

// sort data_list[] into ascendent order
template<typename dataType>
void insertionSortAsc(dataType *data_list, size_t num)
{
	unsigned int i, j;
	dataType tmp;
	for (i = 1; i < num; i++)
		for (j = i; j > 0 && data_list[j - 1] > data_list[j]; j--)
		{
			tmp = data_list[j - 1];
			data_list[j - 1] = data_list[j];
			data_list[j] = tmp;
		}
}

#endif