//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapBlock - Data structures holding
// the result of the alignment of a sequence read
// to a BWT
// 
#include "OverlapBlock.h"
#include "BWTAlgorithms.h"

// Return a pointer to the BWT that should be used to extend the block
// this is the opposite BWT that was used in the backwards search
const BWT* OverlapBlock::getExtensionBWT(const BWT* pBWT, const BWT* pRevBWT) const
{
	if(!flags.isTargetRev())
		return pRevBWT;
	else
		return pBWT;
}

// 
AlphaCount OverlapBlock::getCanonicalExtCount(const BWT* pBWT, const BWT* pRevBWT) const
{
	AlphaCount out = BWTAlgorithms::getExtCount(ranges.interval[1], getExtensionBWT(pBWT, pRevBWT));
	if(flags.isQueryComp())
		out.complement();
	return out;
}


//
void printList(OverlapBlockList* pList)
{
	for(OverlapBlockList::iterator i = pList->begin(); i != pList->end(); ++i)
	{
		std::cout << "Block: " << *i << "\n";
	}
}

//
void removeSubMaximalBlocks(OverlapBlockList* pList)
{
	// This algorithm removes any sub-maximal OverlapBlocks from pList
	// The list is sorted by the left coordinate and iterated through
	// if two adjacent blocks overlap they are split into maximal contiguous regions
	// with resolveOverlap. The resulting list is merged back into pList. This process
	// is repeated until each block in pList is a unique range
	// The bookkeeping in the intersecting case could be more efficient 
	// but the vast vast majority of the cases will not have overlapping 
	// blocks and will just iterate straight through the list so we do
	// something simple here.
	pList->sort(OverlapBlock::sortIntervalLeft);
	OverlapBlockList::iterator iter = pList->begin();
	OverlapBlockList::iterator last = pList->end();
	--last;

	while(iter != pList->end())
	{
		OverlapBlockList::iterator next = iter;
		++next;

		if(next == pList->end())
			break;

		// Check if iter and next overlap
		if(Interval::isIntersecting(iter->ranges.interval[0].lower, iter->ranges.interval[0].upper, 
			                        next->ranges.interval[0].lower, next->ranges.interval[0].upper))
		{
			OverlapBlockList resolvedList = resolveOverlap(*iter, *next);
			
			// Merge the new elements in and start back from the beginning of the list
			pList->erase(iter);
			pList->erase(next);
			pList->merge(resolvedList, OverlapBlock::sortIntervalLeft);
			iter = pList->begin();

			//std::cout << "After splice: \n";
			//printList(pList);
		}
		else
		{
			++iter;
		}
	}
}

// In rare cases, the overlap blocks may represent sub-maximal overlaps between reads
// we need to distinguish these cases and remove the sub-optimal hits. This
// function splits two overlapping OverlapBlocks into up to three distinct
// blocks, keeping the maximal (longest) overlap at each stage.
OverlapBlockList resolveOverlap(const OverlapBlock& A, const OverlapBlock& B)
{
	OverlapBlockList outList;

	// Check if A and B have the same overlap length, if so they must be 
	// identical blocks (resulting from different seeds) and we can remove one
	if(A.overlapLen == B.overlapLen)
	{
		if(A.ranges.interval[0].lower == B.ranges.interval[0].lower &&
		   A.ranges.interval[0].upper == B.ranges.interval[0].upper)
		{
			outList.push_back(A);
			return outList;
		}
		else
		{
			std::cerr << "Error in resolveOverlap: Overlap blocks with same length do not \
			have same coordinates\n";
			assert(false);
		}	
	}
	// A and B must have different overlap lengths
	assert(A.overlapLen != B.overlapLen);

	// Determine which of A and B have a higher overlap
	const OverlapBlock* pHigher;
	const OverlapBlock* pLower;
	if(A.overlapLen > B.overlapLen)
	{
		pHigher = &A;
		pLower = &B;
	}
	else
	{
		pHigher = &B;
		pLower = &A;
	}

	// Complicated logic follows
	// We always want the entirity of the block with the longer
	// overlap so it is added to outList unmodified
	outList.push_back(*pHigher);

	// The lower block can be split into up to two pieces:
	// Case 1:
	//     Lower  ------ 
	//     Higher    ------
	//     Result ---
	//
	// Case 2:
	//     Lower  -----------
	//     Higher    ------
	//     Result ---      --
	//
	// Case 3:
	//     Lower  ------
	//     Higher ------
	//     Result (empty set)

	// It is unclear whether case 2 can happen in reality but we handle it 
	// here anyway. Further complicating matters is that the BWTIntervalPair
	// keeps track of both the BWT coordinates for the backwards search
	// and forward search and we must take care to ensure that both intervals
	// are updated and the mapping between them is sane. 
	//
	// Since the ordering of reads within the two intervals
	// is equal, by symmetrically performing the same operation on both intervals 
	// we preserve the mapping from the forward interval to the reverse. For instance 
	// if we contract the forward interval from [0 2] to [0 1] we must also contract the
	// reverse interval from [37 39] to [37 38].

	// Potentially split the lower block into two blocks
	OverlapBlock split = *pLower;

	// Left-hand split
	if(pLower->ranges.interval[0].lower < pHigher->ranges.interval[0].lower)
	{
		split.ranges.interval[0].lower = pLower->ranges.interval[0].lower;
		split.ranges.interval[0].upper = pHigher->ranges.interval[0].lower - 1; // inclusive
		
		// Calculate the new size of the interval and apply it to the reverse interval
		int64_t diff = split.ranges.interval[0].upper - split.ranges.interval[0].lower;

		split.ranges.interval[1].lower = pLower->ranges.interval[1].lower;
		split.ranges.interval[1].upper = split.ranges.interval[1].lower + diff;

		assert(split.ranges.interval[0].size() == split.ranges.interval[1].size());
		assert(split.ranges.interval[0].isValid());
		assert(split.ranges.interval[1].isValid());
		outList.push_back(split);
	}

	// Right-hand split
	if(pLower->ranges.interval[0].upper > pHigher->ranges.interval[0].upper)
	{
		split.ranges.interval[0].lower = pHigher->ranges.interval[0].upper + 1; // inclusive
		split.ranges.interval[0].upper = pLower->ranges.interval[0].upper;
		
		// Calculate the new size of the interval and apply it to the reverse interval
		int64_t diff = split.ranges.interval[0].upper - split.ranges.interval[0].lower;

		split.ranges.interval[1].upper = pLower->ranges.interval[1].upper;
		split.ranges.interval[1].lower = split.ranges.interval[1].upper - diff; 

		assert(split.ranges.interval[0].size() == split.ranges.interval[1].size());
		assert(split.ranges.interval[0].isValid());
		assert(split.ranges.interval[1].isValid());
		outList.push_back(split);
	}

	if(outList.size() == 3)
	{
		WARN_ONCE("Overlap block was split into 3 segments");
	}

	// Sort the outlist by left coordinate
	outList.sort(OverlapBlock::sortIntervalLeft);
	return outList;
}

