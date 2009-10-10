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

template<typename IterType, typename BucketFunctor>
void histogramSort(IterType start, IterType end, BucketFunctor func, int depth = 0)
{
	// Typedef the base value of the iterator
	typedef typename std::iterator_traits<IterType>::value_type base_value;
	typedef std::list<base_value> base_list;
	typedef typename base_list::iterator list_iterator;

	// Get the number of buckets needed
	int numBuckets = func.getNumBuckets();

	// Set the functor's offset
	func.setBucketDepth(depth);

	// Allocate an array for the bucket start points and the bucket end points
	// TODO: Could remove at least one of these arrays
	size_t* bucket_starts = new size_t[numBuckets];
	size_t* bucket_ends = new size_t[numBuckets];
	size_t* bucket_counts = new size_t[numBuckets];

	for(int i = 0; i < numBuckets; ++i)
	{
		bucket_starts[i] = 0;
		bucket_ends[i] = 0;
		bucket_counts[i] = 0;
	}	

	// Place the items into buckets
	for(IterType iter = start; iter != end; ++iter)
	{
		int bucket_id = func(*iter);
		assert(bucket_id < numBuckets);
		bucket_counts[bucket_id]++;
	}

	// Set up the start/end points
	size_t sum = 0;
	for(int i = 0; i < numBuckets; ++i)
	{
		size_t val = bucket_counts[i];
		bucket_starts[i] = sum;
		sum += val;
		bucket_ends[i] = sum;
		//std::cout << "Bucket count: " << bucket_counts[i] << " start: " << bucket_starts[i] << " end: " << bucket_ends[i] << "\n";
		
	}

	// Now, consider the bucket_starts positions to be the next (potentially) unsorted position
	// for each bucket. Iterate through these values cycles the elements into place
	int curr_bucket = 0; // Start on the first bucket that has element
	while(bucket_starts[curr_bucket] == bucket_ends[curr_bucket])
		++curr_bucket;

	while(curr_bucket != numBuckets)
	{
		//std::cout << "Curr Bucket: " << curr_bucket << "\n";
		//std::cout << "CurrIDX: " << bucket_starts[curr_bucket] << "\n";
		
		int elem_offset = bucket_starts[curr_bucket];
		IterType currIter = start + elem_offset;
		
		// Get the bucket for this element
		int bucket_id = func(*currIter);

		// if the element is in the correct bucket, just advance the pointer
		if(bucket_id == curr_bucket)
			bucket_starts[curr_bucket]++;
		else
		{
			// Cycle the elements until we land back in this bucket, at which point
			// we have sorted one position further in the current bucket

			// Hold the info for the element that has been displaced
			base_value displaced = *currIter;
			bool stop = false;

			while(!stop)
			{
				// Get the offset of the position this element should be in
				int displaced_bucket = func(displaced);

				if(displaced_bucket == curr_bucket)
				{
					// The displaced element belongs back in the original bucket, place it there
					*currIter = displaced;
					bucket_starts[curr_bucket]++;
					stop = true;
				}
				else
				{
					int cycle_offset = bucket_starts[displaced_bucket];
					IterType cycleIter = start + cycle_offset;

					// Make space for the incoming element
					base_value swap = *cycleIter;
					*cycleIter = displaced;
					displaced = swap;
					bucket_starts[displaced_bucket]++;
				}
			}
		}

		// At this point, the bucket_start for the current bucket has moved 1 element
		// Other buckets may have advanced as well
		//
		while(curr_bucket < numBuckets && bucket_starts[curr_bucket] == bucket_ends[curr_bucket])
			++curr_bucket;
	}

	// Recompute the start/end points
	sum = 0;
	for(int i = 0; i < numBuckets; ++i)
	{
		size_t val = bucket_counts[i];
		bucket_starts[i] = sum;
		sum += val;
		bucket_ends[i] = sum;
	}

#if 0
	std::cerr << "Warning: histogram sort validation is turned on\n";
	// Iterate through each bucket, ensuring that each element belongs there
	for(int i = 0; i < numBuckets; ++i)
	{
		for(size_t j = bucket_starts[i]; j != bucket_ends[i]; ++j)
		{
			IterType validation_iter = start + j;
			if(func(*validation_iter) != i)
			{
				std::cerr << "Element " << *validation_iter << " belongs in bucket " << 
					func(*validation_iter) << " not bucket " << i << "\n";
				assert(false);
			}
		}
	}
#endif

	// Finally, sort each bucket
	for(int i = 0; i < numBuckets; ++i)
	{
		IterType bucket_start_iter = start + bucket_starts[i];
		IterType bucket_end_iter = start + bucket_ends[i];
		//std::cout << "Sorting bucket " << i << "\n";
		const int max_depth = 100;
		const size_t min_elements = 100;

		// terminate the bucket sort under 3 conditions:
		//  1) the bucket is degenerate (no more histogram sorting can be done) - this comes from the functor
		//  2) the number of elements in the bucket is low
		//  3) the max depth has been hit
		if(bucket_counts[i] < min_elements || depth >= max_depth || func.isBucketDegenerate(i))
		{
			std::sort(bucket_start_iter, bucket_end_iter, func);
		}
		else
		{
			histogramSort(bucket_start_iter, bucket_end_iter, func, depth + 1);
		}
	}

	// delete arrays
	delete [] bucket_counts;
	delete [] bucket_starts;
	delete [] bucket_ends;
}



#endif
