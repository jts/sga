//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// bucketSort - perform a bucket sort over a range of values
//

#ifndef BUCKETSORT_H
#define BUCKETSORT_H
#include <list>

template<typename IterType, typename BucketFunctor>
void _insertionSort(IterType start, IterType end, BucketFunctor func)
{
	if(start == end)
		return;
	typedef typename std::iterator_traits<IterType>::value_type base_value;
	IterType i = start + 1;
	for(IterType i = start + 1; i != end; ++i)
	{
		base_value v = *i;
		IterType j;
		for(j = i - 1; j >= start; --j)
		{
			bool cmp = !func(v, *j);
			if(cmp)
				break;
			*(j+1) = *j;
		}
		*(j+1) = v;
	}
}

template<typename IterType, typename BucketFunctor>
void bucketSort(IterType start, IterType end, BucketFunctor func)
{

	// Typedef the base value of the iterator
	typedef typename std::iterator_traits<IterType>::value_type base_value;
	typedef std::list<base_value> base_list;
	typedef typename base_list::iterator list_iterator;

	// Get the number of buckets needed
	int numBuckets = func.getNumBuckets();

	// Allocate an array of lists
	base_list* buckets = new base_list[numBuckets];

	// Place the items into buckets
	for(IterType iter = start; iter != end; ++iter)
	{
		int bucket_id = func(*iter);
		assert(bucket_id < numBuckets);
		buckets[bucket_id].push_back(*iter);
	}

	// Place the elements back into the array
	IterType curr = start;
	for(int i = 0; i < numBuckets; ++i)
	{
		IterType bucket_start = curr;
		size_t count = 0;
		for(list_iterator bucket_iter = buckets[i].begin(); bucket_iter != buckets[i].end(); ++bucket_iter)
		{
			++count;
			*curr++ = *bucket_iter;
		}
		std::sort(bucket_start, curr, func);
	}
	assert(curr == end);

	// delete lists
	delete [] buckets;
}

#endif
