//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "OverlapAlgorithm.h"
#include "ASQG.h"
#include <math.h>

// Collect the complete set of overlaps in pOBOut
static const AlignFlags sufPreAF(false, false, false);
static const AlignFlags prePreAF(false, true, true);
static const AlignFlags sufSufAF(true, false, true);
static const AlignFlags preSufAF(true, true, false);

//#define TEMPDEBUG 1

// Perform the overlap
OverlapResult OverlapAlgorithm::overlapRead(const SeqRecord& read, OverlapBlockList* pOutList) const
{
	OverlapResult r;
	r = overlapReadInexact(read, pOutList);
	//r = overlapReadExact(read, pOutList);
	
	// If the read is a substring of some other read, clear its overlap block list
	if(r.isSubstring)
		pOutList->clear();
	return r;
}

//
OverlapResult OverlapAlgorithm::overlapReadInexact(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();

#ifdef TEMPDEBUG
	std::cout << "Overlapping read " << read.id << " suffix\n";
#endif

	// Match the suffix of seq to prefixes
	findOverlapBlocksInexact(seq, m_pBWT, m_pRevBWT, sufPreAF, &obWorkingList, pOBOut, result);
	findOverlapBlocksInexact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, &obWorkingList, pOBOut, result);

	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
		obWorkingList.clear();
	}
	else
	{
		pOBOut->splice(pOBOut->end(), obWorkingList);
		assert(obWorkingList.empty());
	}

#ifdef TEMPDEBUG
	std::cout << "Overlapping read " << read.id << " prefix\n";
#endif

	// Match the prefix of seq to suffixes
	findOverlapBlocksInexact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, &obWorkingList, pOBOut, result);
	findOverlapBlocksInexact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, &obWorkingList, pOBOut, result);

	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
		obWorkingList.clear();
	}
	else
	{
		pOBOut->splice(pOBOut->end(), obWorkingList);
		assert(obWorkingList.empty());
	}

	return result;
}

// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
OverlapResult OverlapAlgorithm::overlapReadExact(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;

	// The complete set of overlap blocks are collected in obWorkingList
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obWorkingList;
	std::string seq = read.seq.toString();

	// Irreducible overlaps only
	WARN_ONCE("Irreducible-only assumptions: All reads are the same length")

	// Match the suffix of seq to prefixes
	findOverlapBlocksExact(seq, m_pBWT, m_pRevBWT, sufPreAF, &obWorkingList, pOBOut, result);
	findOverlapBlocksExact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, &obWorkingList, pOBOut, result);

	//	
	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
		obWorkingList.clear();
	}
	else
	{
		pOBOut->splice(pOBOut->end(), obWorkingList);
		assert(obWorkingList.empty());
	}
	// Match the prefix of seq to suffixes
	findOverlapBlocksExact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, &obWorkingList, pOBOut, result);
	findOverlapBlocksExact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, &obWorkingList, pOBOut, result);

	//
	if(m_bIrreducible)
	{
		computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
		obWorkingList.clear();
	}
	else
	{
		pOBOut->splice(pOBOut->end(), obWorkingList);
		assert(obWorkingList.empty());
	}

	return result;
}

// Write overlap results to an ASQG file
void OverlapAlgorithm::writeResultASQG(std::ostream& writer, const SeqRecord& read, const OverlapResult& result) const
{
	ASQG::VertexRecord record(read.id, read.seq.toString());
	record.setSubstringTag(result.isSubstring);
	record.write(writer);
}

// Write overlap blocks out to a file
void OverlapAlgorithm::writeOverlapBlocks(std::ostream& writer, size_t readIdx, const OverlapBlockList* pList) const
{
	// Write the hits to the file
	if(!pList->empty())
	{
		// Write the header info
		size_t numBlocks = pList->size();
		writer << readIdx << " " << numBlocks << " ";
		//std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
		for(OverlapBlockList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
		{
			writer << *iter << " ";
		}
		writer << "\n";
	}
}

// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
void OverlapAlgorithm::findOverlapBlocksExact(const std::string& w, const BWT* pBWT,
      									      const BWT* pRevBWT, const AlignFlags& af, 
											  OverlapBlockList* pOBList, OverlapBlockList* pOBFinal, 
											  OverlapResult& result) const
{
	// All overlaps are added to this list and then sub-maximal overlaps are removed
	OverlapBlockList workingList;

	// The algorithm is as follows:
	// We perform a backwards search using the FM-index for the string w.
	// As we perform the search we collect the intervals 
	// of the significant prefixes (len >= minOverlap) that overlap w.
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
	
	// Collect the OverlapBlocks
	for(size_t i = start - 1; i >= 1; --i)
	{
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
		size_t overlapLen = l - i;
		if(overlapLen >= m_minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal
			// These are the proper prefixes (they are the start of a read)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// The probe interval contains the range of proper prefixes
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				workingList.push_back(OverlapBlock(probe, overlapLen, 0, af));
			}
		}
	}

	// Determine if this sequence is contained and should not be processed further
	BWTAlgorithms::updateBothL(ranges, w[0], pBWT);

	// Ranges now holds the interval for the full-length read
	// To handle containments, we output the overlapBlock to the final overlap block list
	// and it will be processed later
	// Two possible containment cases:
	// 1) This read is a substring of some other read
	// 2) This read is identical to some other read
	
	// Case 1 is indicated by the existance of a non-$ left or right hand extension
	// In this case we return no alignments for the string
	AlphaCount left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], pBWT);
	AlphaCount right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);
	if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
	{
		result.isSubstring = true;
	}
	else
	{
		BWTAlgorithms::updateBothL(ranges, '$', pBWT);
		if(ranges.isValid())
			pOBFinal->push_back(OverlapBlock(ranges, w.length(), 0, af));
	}

	// Remove sub-maximal OverlapBlocks and move the remainder to the output list
	removeSubMaximalBlocks(&workingList);
	pOBList->splice(pOBList->end(), workingList);
	return;
}

// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive
void OverlapAlgorithm::findOverlapBlocksInexact(const std::string& w, const BWT* pBWT, 
                                                const BWT* pRevBWT, const AlignFlags& af, 
												OverlapBlockList* pOBList, OverlapBlockList* pOBFinal, 
												OverlapResult& result) const
{
	int len = w.length();
	int overlap_region_left = len - m_minOverlap;
	SearchSeedVector* pCurrVector = new SearchSeedVector;
	SearchSeedVector* pNextVector = new SearchSeedVector;
	OverlapBlockList partialWorkingList;
	OverlapBlockList fullWorkingList;
	SearchSeedVector::iterator iter;

	// Create and extend the initial seeds
	int actual_seed_length = m_seedLength;
	int actual_seed_stride = m_seedStride;
	if(actual_seed_length == 0)
	{
		// Calculate a seed length and stride that will guarantee all overlaps
		// with error rate m_errorRate will be found
		calculateSeedParameters(w, actual_seed_length, actual_seed_stride);
	}

	createSearchSeeds(w, pBWT, pRevBWT, actual_seed_length, actual_seed_stride, pCurrVector);
	extendSeedsExactRight(w, pBWT, pRevBWT, ED_RIGHT, pCurrVector, pNextVector);
	pCurrVector->clear();
	pCurrVector->swap(*pNextVector);
	assert(pNextVector->empty());

	int num_steps = 0;

	// Perform the inexact extensions
	while(!pCurrVector->empty())
	{
		iter = pCurrVector->begin();
		while(iter != pCurrVector->end())
		{
			SearchSeed& align = *iter;

			// If the current aligned region is right-terminal
			// and the overlap is greater than minOverlap, try to find overlaps
			// or containments
			if(align.right_index == len - 1)
			{
				double align_error = align.calcErrorRate();

				// Check for overlaps
				if(align.left_index <= overlap_region_left && align_error <= m_errorRate)
				{
					int overlapLen = len - align.left_index;
					BWTIntervalPair probe = align.ranges;
					BWTAlgorithms::updateBothL(probe, '$', pBWT);
					
					// The probe interval contains the range of proper prefixes
					if(probe.interval[1].isValid())
					{
						assert(probe.interval[1].lower > 0);
						OverlapBlock nBlock(OverlapBlock(probe, overlapLen, align.z, af, align.history));
						if(overlapLen == len)
							fullWorkingList.push_back(nBlock);
						else
							partialWorkingList.push_back(nBlock);
					}
				}

				// Check for containments
				// If the seed is left-terminal and there are [ACGT] left/right extensions of the sequence
				// this read must be a substring of another read
				if(align.left_index == 0)
				{
					AlphaCount left_ext = BWTAlgorithms::getExtCount(align.ranges.interval[0], pBWT);
					AlphaCount right_ext = BWTAlgorithms::getExtCount(align.ranges.interval[1], pRevBWT);
					if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
						result.isSubstring = true;
				}
			}

			// Extend the seed to the right/left
			if(align.dir == ED_RIGHT)
				extendSeedInexactRight(align, w, pBWT, pRevBWT, pNextVector);
			else
                extendSeedInexactLeft(align, w, pBWT, pRevBWT, pNextVector);
			++iter;
			//pCurrVector->erase(iter++);
		}
		pCurrVector->clear();
		assert(pCurrVector->empty());
		pCurrVector->swap(*pNextVector);

		// Remove identical seeds after we have performed seed_len steps
		// as there is now the chance of identical seeds
		if(num_steps % actual_seed_stride == 0)
		{
			std::sort(pCurrVector->begin(), pCurrVector->end(), SearchSeed::compareLeftRange);
			SearchSeedVector::iterator end_iter = std::unique(pCurrVector->begin(), pCurrVector->end(), 
			                                                       SearchSeed::equalLeftRange);
			pCurrVector->resize(end_iter - pCurrVector->begin());
		}
		++num_steps;
	}

	WARN_ONCE("TODO: Collect all blocks in the same list then split AFTER removing submaximal");
	// parse the full working list, which has containment overlaps
	removeSubMaximalBlocks(&fullWorkingList);
	pOBFinal->splice(pOBFinal->end(), fullWorkingList);

	// parse the partial block working list, which has the proper overlaps
	removeSubMaximalBlocks(&partialWorkingList);
	pOBList->splice(pOBList->end(), partialWorkingList);

	delete pCurrVector;
	delete pNextVector;
}

// Calculate the seed length and stride to ensure that we will find all 
// all overlaps with error rate at most m_errorRate.
// To overlap two reads allowing for d errors, we create d+1 seeds so that at least one seed will match
// exactly. d is a function of the overlap length so we define the seed length using the minimum overlap
// parameter. We then tile seeds across the read starting from the end such that for every overlap length
// x, there are at least floor(x * error_rate) + 1 seeds.
void OverlapAlgorithm::calculateSeedParameters(const std::string& w, int& seed_length, int& seed_stride) const
{
	int read_len = w.length();
	seed_length = 0;
	
	// The maximum possible number of differences occurs for a fully-aligned read
	int max_diff_high = static_cast<int>(m_errorRate * read_len);

	// Calculate the seed length to use
	// If the error rate is so low that no differences are possible just seed
	// over the entire minOverlap region
	if(max_diff_high > 0)
	{
		// Calculate the maximum number of differences between two sequences that overlap
		// by minOverlap
		int max_diff_low = static_cast<int>(m_errorRate * m_minOverlap);

		// Calculate the length of the region such that max_diff_low / region <= error_rate
		// If max_diff_low is zero, calculate the minimal region to find one error
		if(max_diff_low == 0)
			max_diff_low = 1;
		
		int seed_region_length = ceil(max_diff_low / m_errorRate);
		int num_seeds_low = max_diff_low + 1;
		seed_length = static_cast<int>(seed_region_length / num_seeds_low);
		assert(seed_length <= static_cast<int>(m_minOverlap));
	}
	else
	{
		seed_length = m_minOverlap;
	}
	seed_stride = seed_length;	

}

// Create and intialize the search seeds
int OverlapAlgorithm::createSearchSeeds(const std::string& w, const BWT* pBWT, 
                                        const BWT* pRevBWT, int seed_length, int seed_stride,
										SearchSeedVector* pOutVector) const
{
	// The maximum possible number of differences occurs for a fully-aligned read
	int read_len = w.length();
	int max_diff_high = static_cast<int>(m_errorRate * read_len);

	static int once = 1;
	if(once)
	{
		printf("Using seed length %d, seed stride %d, max diff %d\n", seed_length, seed_stride, max_diff_high);
		once = 0;
	}	

	// Start the seeds at the end of the read
	int seed_start = read_len - seed_length;

	while(seed_start >= 0)
	{
		SearchSeed seed;
		seed.left_index = seed_start;
		seed.right_index = seed_start;
		seed.dir = ED_RIGHT;
		seed.seed_len = seed_length;
		seed.z = 0;
		seed.maxDiff = max_diff_high;

		// Initialize the left and right suffix array intervals
		char b = w[seed.left_index];
		BWTAlgorithms::initIntervalPair(seed.ranges, b, pBWT, pRevBWT);		
		pOutVector->push_back(seed);

		seed_start -= seed_stride;

		// Only create one seed in the exact case
		if(max_diff_high == 0)
			break;
	}
	return seed_length;
}

// Extend all the seeds in pInVector to the right over the entire seed range
void OverlapAlgorithm::extendSeedsExactRight(const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT,
                                             ExtendDirection /*dir*/, const SearchSeedVector* pInVector, 
											 SearchSeedVector* pOutVector) const
{
	for(SearchSeedVector::const_iterator iter = pInVector->begin(); iter != pInVector->end(); ++iter)
	{
		SearchSeed align = *iter;
		bool valid = true;
		while(align.isSeed())
		{
			++align.right_index;
			char b = w[align.right_index];
			BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
			if(!align.isIntervalValid(RIGHT_INT_IDX))
			{
				valid = false;
				break;
			}
		}
		
		if(valid)
			pOutVector->push_back(align);
	}
}

//
void OverlapAlgorithm::extendSeedInexactRight(SearchSeed& seed, const std::string& w, const BWT* /*pBWT*/, 
                                              const BWT* pRevBWT, SearchSeedVector* pOutVector) const
{
	// If this alignment has run all the way to the end of the sequence
	// switch it to be a left extension sequence
	if(seed.right_index == int(w.length() - 1))
	{
		seed.dir = ED_LEFT;
		pOutVector->push_back(seed);
		return;
	}

	++seed.right_index;
	
	if(!seed.allowMismatches())
	{
		char b = w[seed.right_index];
		BWTAlgorithms::updateBothR(seed.ranges, b, pRevBWT);
		if(seed.isIntervalValid(RIGHT_INT_IDX))
			pOutVector->push_back(seed);
	}
	else
	{
		for(int i = 0; i < 4; ++i)
		{
			char b = ALPHABET[i];
			SearchSeed branched = seed;
			BWTAlgorithms::updateBothR(branched.ranges, b, pRevBWT);
			if(branched.isIntervalValid(RIGHT_INT_IDX))
			{
				if(b != w[seed.right_index])
				{
					++branched.z;
					// The history coordinates are wrt the right end of the read
					// so that each position corresponds to the length of the overlap
					// including that position
					branched.history.add(w.length() - seed.right_index, b);
				}
				pOutVector->push_back(branched);
			}
		}
	}
}

//
void OverlapAlgorithm::extendSeedInexactLeft(SearchSeed& seed, const std::string& w, 
                                             const BWT* pBWT, const BWT* /*pRevBWT*/,
                                             SearchSeedVector* pOutVector) const
{
	//printf("ProcessingLEFT: "); align.print(w);
	--seed.left_index;
	if(seed.left_index >= 0)
	{
		if(!seed.allowMismatches())
		{
			// Extend exact
			char b = w[seed.left_index];
			BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
			if(seed.isIntervalValid(LEFT_INT_IDX))
				pOutVector->push_back(seed);
		}
		else
		{
			for(int i = 0; i < 4; ++i)
			{
				char b = ALPHABET[i];
				SearchSeed branched = seed;
				BWTAlgorithms::updateBothL(branched.ranges, b, pBWT);
				if(branched.isIntervalValid(LEFT_INT_IDX))
				{
					if(b != w[seed.left_index])
					{
						++branched.z;
						// The history coordinates are wrt the right end of the read
						// so that each position corresponds to the length of the overlap
						// including that position						
						branched.history.add(w.length() - seed.left_index, b);
					}
					pOutVector->push_back(branched);
				}
			}
		}
	}
}


// Calculate the irreducible blocks from the vector of OverlapBlocks
void OverlapAlgorithm::computeIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, 
                                                OverlapBlockList* pOBList, 
												OverlapBlockList* pOBFinal) const
{
	// processIrreducibleBlocks requires the pOBList to be sorted in descending order
	pOBList->sort(OverlapBlock::sortSizeDescending);
	_processIrreducibleBlocksInexact(pBWT, pRevBWT, *pOBList, pOBFinal);
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::_processIrreducibleBlocksExact(const BWT* pBWT, const BWT* pRevBWT, 
                                                      OverlapBlockList& obList, 
													  OverlapBlockList* pOBFinal) const
{
	if(obList.empty())
		return;
	
	// Count the extensions in the top level (longest) blocks first
	int topLen = obList.front().overlapLen;
	AlphaCount ext_count;
	OBLIter iter = obList.begin();
	while(iter != obList.end() && iter->overlapLen == topLen)
	{
		ext_count += iter->getCanonicalExtCount(pBWT, pRevBWT);
		++iter;
	}
	
	// Three cases:
	// 1) The top level block has ended as it contains the extension $. Output TLB and end.
	// 2) There is a singular unique extension base for all the blocks. Update all blocks and recurse.
	// 3) There are multiple extension bases, branch and recurse.
	// If some block other than the TLB ended, it must be contained within the TLB and it is not output
	// or considered further. 
	// Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
	// must have the same sequence as the top one and are hence considered to be contained with the top element.
	if(ext_count.get('$') > 0)
	{
		// An irreducible overlap has been found. It is possible that there are two top level blocks
		// (one in the forward and reverse direction). Since we can't decide which one
		// contains the other at this point, we output hits to both. Under a fixed 
		// length string assumption one will be contained within the other and removed later.
		OBLIter tlbIter = obList.begin();
		while(tlbIter != obList.end() && tlbIter->overlapLen == topLen)
		{
			// Ensure the tlb is actually terminal and not a substring block
			AlphaCount test_count = tlbIter->getCanonicalExtCount(pBWT, pRevBWT);
			assert(test_count.get('$') > 0);
			pOBFinal->push_back(OverlapBlock(*tlbIter));
			++tlbIter;
		} 
		return;
	}
	else
	{
		// Count the rest of the blocks
		while(iter != obList.end())
		{
			ext_count += iter->getCanonicalExtCount(pBWT, pRevBWT);
			++iter;
		}

		if(ext_count.hasUniqueDNAChar())
		{

			// Update all the blocks using the unique extension character
			// This character is in the canonical representation wrt to the query
			char b = ext_count.getUniqueDNAChar();
			updateOverlapBlockRangesRight(pBWT, pRevBWT, obList, b);
			return _processIrreducibleBlocksExact(pBWT, pRevBWT, obList, pOBFinal);
		}
		else
		{
			for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
			{
				char b = ALPHABET[idx];
				if(ext_count.get(b) > 0)
				{
					OverlapBlockList branched = obList;
					updateOverlapBlockRangesRight(pBWT, pRevBWT, branched, b);
					_processIrreducibleBlocksExact(pBWT, pRevBWT, branched, pOBFinal);
				}
			}
		}
	}
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::_processIrreducibleBlocksInexact(const BWT* pBWT, const BWT* pRevBWT, 
                                                      OverlapBlockList& obList, 
													  OverlapBlockList* pOBFinal) const
{
	OverlapBlockList* pOBNext = new OverlapBlockList;

	if(obList.empty())
		return;

	// Count the extensions in the top level (longest) blocks first
	int extension_length = 0;
	bool done = false;

	while(!done)
	{
		int topLen = obList.front().overlapLen;
		AlphaCount ext_count;
		OBLIter iter = obList.begin();
		while(iter != obList.end() && iter->overlapLen == topLen)
		{
			ext_count += iter->getCanonicalExtCount(pBWT, pRevBWT);
			++iter;
		}
		
		if(ext_count.get('$') > 0)
		{
			//std::cout << "TLB of length " << topLen << " has ended\n";
			// The end of the top level block(s) have been found
			OBLIter tlbIter = obList.begin();
			while(tlbIter != obList.end() && tlbIter->overlapLen == topLen)
			{
				AlphaCount test_count = tlbIter->getCanonicalExtCount(pBWT, pRevBWT);
				assert(test_count.get('$') > 0);
				
				// If this block has been marked eliminated, it is transitive wrt this read
				// and we do not output the block
				if(!tlbIter->isEliminated)
				{
#ifdef TEMPDEBUG
					std::cout << "Pushing block of length " << tlbIter->overlapLen << "\n";
#endif
					pOBFinal->push_back(*tlbIter);
					tlbIter->isEliminated = true;
				}
				
				// Mark any blocks that are shorter than the TLB as transitive
				// if the inferred error rate is less than the threshold
				OBLIter transIter = tlbIter;
				++transIter;
				for(; transIter != obList.end(); ++transIter)
				{
					if(transIter->overlapLen == topLen)
					{
						WARN_ONCE("Skipping identical-length block");
						continue;
					}
					else if(!transIter->isEliminated)
					{
						// Compute error rate between the transIter block and tlbIter block,
						// mark all the blocks that have error rate wrt the tlb lower than 
						// m_errorRate as eliminated as they must be transitive edges
						int backwards_diff = SearchHistory::countDifferences(tlbIter->backHistory, transIter->backHistory, transIter->overlapLen);
						int forward_diff = SearchHistory::countDifferences(tlbIter->forwardHistory, transIter->forwardHistory, extension_length);
						int trans_overlap_length = transIter->overlapLen + extension_length;
						double er = static_cast<double>(backwards_diff + forward_diff) / trans_overlap_length;
						
#ifdef TEMPDEBUG
						std::cout << "OL: " << transIter->overlapLen << "\n";
						std::cout << "TLB BH: " << tlbIter->backHistory << "\n";
						std::cout << "TB  BH: " << transIter->backHistory << "\n";
						std::cout << "TLB FH: " << tlbIter->forwardHistory << "\n";
						std::cout << "TB  FH: " << transIter->forwardHistory << "\n";
						std::cout << "IOL: " << trans_overlap_length << " TD: " << (backwards_diff + forward_diff) << "\n";
						std::cout << "Block of length " << transIter->overlapLen << " has ier: " << er << "\n";
#endif
						// 
						if(er <= m_errorRate)
						{
							//std::cout << "Marking block of length " << transIter->overlapLen << " as eliminated\n";
							transIter->isEliminated = true;
						}
					}

					// Shift the block to the next list
					pOBNext->push_back(*transIter);
				}

				++tlbIter;
			} 
	
			// Check if all the blocks have been eliminated
			bool all_eliminated = true;
			for(OBLIter nextIter = pOBNext->begin(); nextIter != pOBNext->end(); ++nextIter)
			{
				if(!nextIter->isEliminated)
				{
					all_eliminated = false;
					break;
				}
			}
			done = all_eliminated;
		}
		else
		{
			// Extend all the blocks one base
			++extension_length;
			for(OBLIter iter = obList.begin(); iter != obList.end(); ++iter)
			{
				for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
				{
					OverlapBlock branched = *iter;
					char b = ALPHABET[idx];
					char cb = iter->flags.isQueryComp() ? complement(b) : b;
					BWTAlgorithms::updateBothR(branched.ranges, cb, branched.getExtensionBWT(pBWT, pRevBWT));
					branched.forwardHistory.add(extension_length, b);
					if(branched.ranges.isValid())
					{
						pOBNext->push_back(branched);
					}
				}
			}
		}

		// All the remaining overlap blocks have been moved to the pOBNext list
		// Clear the current list and swap 
		obList.clear();
		obList.swap(*pOBNext);
	}
	delete pOBNext;
}


// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void OverlapAlgorithm::updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, 
		                                             OverlapBlockList& obList, char b) const
{
	OverlapBlockList::iterator iter = obList.begin(); 
	while(iter != obList.end())
	{
		char cb = iter->flags.isQueryComp() ? complement(b) : b;
		BWTAlgorithms::updateBothR(iter->ranges, cb, iter->getExtensionBWT(pBWT, pRevBWT));
		// remove the block from the list if its no longer valid
		if(!iter->ranges.isValid())
			iter = obList.erase(iter);
		else
			++iter;
	}
}

