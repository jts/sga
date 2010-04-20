//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "OverlapAlgorithm.h"
#include "ASQG.h"
#include <math.h>

#define SWAP_LIST(x, y) pSwap = (x); (x) = (y); (y) = pSwap;

// Collect the complete set of overlaps in pOBOut
static const AlignFlags sufPreAF(false, false, false);
static const AlignFlags prePreAF(false, true, true);
static const AlignFlags sufSufAF(true, false, true);
static const AlignFlags preSufAF(true, true, false);


// Perform the overlap
OverlapResult OverlapAlgorithm::overlapRead(const SeqRecord& read, OverlapBlockList* pOutList) const
{
	OverlapResult r;
	if(!m_bIrreducible)
	{
		//overlapReadExhaustive(read, pOutList);
		r = overlapReadInexact(read, pOutList);
	}
	else
	{
		r = overlapReadIrreducible(read, pOutList);
	}
	
	// If the read is a substring of some other read, clear its overlap block list
	if(r.isSubstring)
		pOutList->clear();
	return r;
}

//
OverlapResult OverlapAlgorithm::overlapReadExhaustive(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	std::string seq = read.seq.toString();

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, pOBOut, pOBOut, result);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, pOBOut, pOBOut, result);

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, pOBOut, pOBOut, result);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, pOBOut, pOBOut, result);
	return result;
}

OverlapResult OverlapAlgorithm::overlapReadInexact(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	std::string seq = read.seq.toString();

	// Match the suffix of seq to prefixes
	overlapSuffixInexact(seq, m_pBWT, m_pRevBWT, m_errorRate, m_minOverlap, sufPreAF, pOBOut, pOBOut, result);
	overlapSuffixInexact(complement(seq), m_pRevBWT, m_pBWT, m_errorRate, m_minOverlap, prePreAF, pOBOut, pOBOut, result);

	// Match the prefix of seq to suffixes
	overlapSuffixInexact(reverseComplement(seq), m_pBWT, m_pRevBWT, m_errorRate, m_minOverlap, sufSufAF, pOBOut, pOBOut, result);
	overlapSuffixInexact(reverse(seq), m_pRevBWT, m_pBWT, m_errorRate, m_minOverlap, preSufAF, pOBOut, pOBOut, result);

	return result;
}

// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
OverlapResult OverlapAlgorithm::overlapReadIrreducible(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
	OverlapResult result;
	// The complete set of overlap blocks are collected in obTemp
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obTemp;
	std::string seq = read.seq.toString();

	// Irreducible overlaps only
	WARN_ONCE("Irreducible-only assumptions: All reads are the same length")

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, &obTemp, pOBOut, result);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, &obTemp, pOBOut, result);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
	obTemp.clear();

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, &obTemp, pOBOut, result);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, &obTemp, pOBOut, result);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);

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
void OverlapAlgorithm::findOverlapBlocks(const std::string& w, const BWT* pBWT,
									  const BWT* pRevBWT, int minOverlap,
									  const AlignFlags& af, OverlapBlockList* pOBList,
                                      OverlapBlockList* pOBFinal, OverlapResult& result)
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
	for(int i = start - 1; i >= 1; --i)
	{
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
		int overlapLen = l - i;
		if(overlapLen >= minOverlap)
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
		// This case isn't handled yet
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

// Set up the alignment blocks and call the alignment function on each block
int OverlapAlgorithm::overlapSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                           double error_rate, int minOverlap, const AlignFlags& af, 
										   OverlapBlockList* pOBList, OverlapBlockList* pOBFinal, OverlapResult& result)
{
	if((int)w.length() < minOverlap)
	{
		std::cerr << "Warning string " << w << " is shorter than minOverlap, it will not be aligned\n";
		return 0;
	}

	// Calculate the number of blocks and the size of each block
	int len = static_cast<int>(w.length());
	
	// Calculate the shortest overlap s.t. 1 difference in the overlap region is less than error_rate
	int max_differences = (int)floor(error_rate * (double)len);

	// Align the reads in segments of shortest_overlap
	// Each segment allows 1 more error than the previous
	int cost = 0;
	int covered = 0;
	for(int i = 0; i <= max_differences; ++i)
	{
		// Calculate the minimum amount of overlap s.t. i differences is less than error_rate
		int min_overlap_size = std::max((int)ceil(i / error_rate), minOverlap);

		// Calculate the endpoint of the block, it is the amount of overlap s.t. i + 1 is less than the error rate
		int max_overlap_size = std::min((int)ceil((i + 1) / error_rate) - 1, len);
		// If the max overlap in this block is less than the minOverlap parameter, skip the block
		if(max_overlap_size < minOverlap)
			continue;

		assert(double(i) / min_overlap_size <= error_rate);
		assert(double(i) / max_overlap_size <= error_rate);
		assert(double(i+1) / max_overlap_size > error_rate);

		// Convert the overlap sizes to block positions, which start from the end of w
		int block_start = len - max_overlap_size;
		int block_end = len - min_overlap_size;
		covered += (block_end - block_start + 1);
		cost += _alignSegmentSimple(w, block_start, block_end, pBWT, pRevBWT, i, af, pOBList, pOBFinal, result);
	}
	assert(covered == (len - minOverlap) + 1);
	return cost;
}

// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive
int OverlapAlgorithm::_alignSegmentSimple(const std::string& w, int block_start, int block_end,
                                           const BWT* pBWT, const BWT* pRevBWT, int maxDiff, 
										   const AlignFlags& af, OverlapBlockList* pOBList, 
										   OverlapBlockList* pOBFinal, OverlapResult& result)
{
	int len = w.length();
	OverlapSeedList* pCurrList = new OverlapSeedList;
	OverlapSeedList* pNextList = new OverlapSeedList;
	OverlapBlockList partialWorkingList;
	OverlapBlockList fullWorkingList;
	OverlapSeedList::iterator iter;

	// Create and extend the initial seeds
	int seed_len = createOverlapSeeds(w, pBWT, pRevBWT, block_start, block_end, maxDiff, pCurrList);
	extendSeedsExactRight(w, pBWT, pRevBWT, ED_RIGHT, pCurrList, pNextList);
	pCurrList->clear();
	pCurrList->swap(*pNextList);
	assert(pNextList->empty());

	int num_steps = 0;

	// Perform the inexact extensions
	while(!pCurrList->empty())
	{
		iter = pCurrList->begin();
		while(iter != pCurrList->end())
		{
			OverlapSeed& align = *iter;

			// Check for overlaps
			if(align.right_index == len - 1)
			{
				// Output overlaps
				if(align.left_index <= block_end)
				{
					int overlapLen = len - align.left_index;
					BWTIntervalPair probe = align.ranges;
					BWTAlgorithms::updateBothL(probe, '$', pBWT);
					
					// The probe interval contains the range of proper prefixes
					if(probe.interval[1].isValid())
					{
						assert(probe.interval[1].lower > 0);
						OverlapBlock nBlock(OverlapBlock(probe, overlapLen, maxDiff - align.z, af));
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
				extendSeedInexactRight(align, w, pBWT, pRevBWT, pNextList);
			else
                extendSeedInexactLeft(align, w, pBWT, pRevBWT, block_start, pNextList);

			pCurrList->erase(iter++);
		}
		assert(pCurrList->empty());
		pCurrList->swap(*pNextList);

		// Remove identical seeds after we have performed seed_len steps
		// as there is now the chance of identical seeds
		++num_steps;
		if(num_steps <= block_end && num_steps % seed_len == 0)
		{
			pCurrList->sort(OverlapSeed::compareLeftRange);
			pCurrList->unique(OverlapSeed::equalLeftRange);
		}
	}

	// parse the full working list, which has containment overlaps
	removeSubMaximalBlocks(&fullWorkingList);
	pOBFinal->splice(pOBFinal->end(), fullWorkingList);

	// parse the partial block working list, which has the proper overlaps
	removeSubMaximalBlocks(&partialWorkingList);
	pOBList->splice(pOBList->end(), partialWorkingList);

	delete pCurrList;
	delete pNextList;
	return 0;
}

int OverlapAlgorithm::createOverlapSeeds(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                          int /*block_start*/, int block_end, int maxDiff, 
										  OverlapSeedList* pOutList)
{
	int len = w.length();
	int seed_start = block_end;
	int seed_end = len;
	int num_seeds = maxDiff + 1;
	int seed_len = (seed_end - seed_start) / num_seeds;

	//printf("srs: %d sre: %d bs: %d be: %d d: %d\n", seed_start, seed_end, block_start, block_end, maxDiff);
	// Populate the initial seeds
	for(int i = 0; i < num_seeds; ++i)
	{
		int pos = seed_start + i*seed_len;
		OverlapSeed align;
		align.left_index = pos;
		align.right_index = pos;
		align.dir = ED_RIGHT;
		align.z = maxDiff;
		align.seed_len = std::min(seed_len, len - pos);
		assert(pos < len);

		//printf("Creating seed at %d to %d, str: %s\n", align.left_index, align.left_index + align.seed_len, w.substr(pos, align.seed_len).c_str());

		// Initialize the left and right suffix array intervals
		char b = w[pos];
		BWTAlgorithms::initIntervalPair(align.ranges, b, pBWT, pRevBWT);		
		pOutList->push_back(align);
	}
	return seed_len;
}

//
void OverlapAlgorithm::extendSeedsExactRight(const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT,
                                             ExtendDirection /*dir*/, const OverlapSeedList* pInList, 
											 OverlapSeedList* pOutList)
{
	for(OverlapSeedList::const_iterator iter = pInList->begin(); iter != pInList->end(); ++iter)
	{
		OverlapSeed align = *iter;
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
			pOutList->push_back(align);
	}
}

//
void OverlapAlgorithm::extendSeedInexactRight(OverlapSeed& seed, const std::string& w, const BWT* /*pBWT*/, 
                                              const BWT* pRevBWT, OverlapSeedList* pOutList)
{
	// If this alignment has run all the way to the end of the sequence
	// switch it to be a left extension sequence
	if(seed.right_index == int(w.length() - 1))
	{
		seed.dir = ED_LEFT;
		pOutList->push_back(seed);
		return;
	}

	++seed.right_index;
	
	if(seed.z == 0)
	{
		char b = w[seed.right_index];
		BWTAlgorithms::updateBothR(seed.ranges, b, pRevBWT);
		if(seed.isIntervalValid(RIGHT_INT_IDX))
			pOutList->push_back(seed);
	}
	else
	{
		for(int i = 0; i < 4; ++i)
		{
			char b = ALPHABET[i];
			OverlapSeed branched = seed;
			BWTAlgorithms::updateBothR(branched.ranges, b, pRevBWT);
			if(branched.isIntervalValid(RIGHT_INT_IDX))
			{
				if(b != w[seed.right_index])
					--branched.z;
				pOutList->push_back(branched);
			}
		}
	}
}

void OverlapAlgorithm::extendSeedInexactLeft(OverlapSeed& seed, const std::string& w, 
                                             const BWT* pBWT, const BWT* /*pRevBWT*/,
                                             int block_start, OverlapSeedList* pOutList)
{
	//printf("ProcessingLEFT: "); align.print(w);
	--seed.left_index;
	if(seed.left_index >= block_start)
	{
		if(seed.z == 0)
		{
			// Extend exact
			char b = w[seed.left_index];
			BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
			if(seed.isIntervalValid(LEFT_INT_IDX))
				pOutList->push_back(seed);
		}
		else
		{
			for(int i = 0; i < 4; ++i)
			{
				char b = ALPHABET[i];
				OverlapSeed branched = seed;
				BWTAlgorithms::updateBothL(branched.ranges, b, pBWT);
				if(branched.isIntervalValid(LEFT_INT_IDX))
				{
					if(b != w[seed.left_index])
						--branched.z;
					pOutList->push_back(branched);
				}
			}
		}
	}
}


// Calculate the irreducible hits from the vector of OverlapBlocks
void OverlapAlgorithm::calculateIrreducibleHits(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList* pOBList, OverlapBlockList* pOBFinal)
{
	// processIrreducibleBlocks requires the pOBList to be sorted in descending order
	pOBList->sort(OverlapBlock::sortSizeDescending);
	processIrreducibleBlocks(pBWT, pRevBWT, *pOBList, pOBFinal);
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::processIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, OverlapBlockList* pOBFinal)
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
	// 1) The top level block has ended, it contains the extension $. Output TLB and end.
	// 2) There is a singular unique extension base for all the blocks. Update all blocks and recurse.
	// 3) There are multiple extension bases, branch and recurse.
	// If some block other than the TLB ended, it must be contained within the TLB and it is not output
	// or considered further. Containments are handled elsewhere.
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
			return processIrreducibleBlocks(pBWT, pRevBWT, obList, pOBFinal);
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
					processIrreducibleBlocks(pBWT, pRevBWT, branched, pOBFinal);
				}
			}
		}
	}
}

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void OverlapAlgorithm::updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, char b)
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

