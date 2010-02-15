//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "OverlapAlgorithm.h"

// Perform the overlap
void OverlapAlgorithm::overlapRead(const SeqItem& read, OverlapBlockList* pOutList) const
{
	if(!m_bIrreducible)
	{
		overlapReadExhaustive(read, pOutList);
	}
	else
	{
		 overlapReadIrreducible(read, pOutList);
	}
}

//
void OverlapAlgorithm::overlapReadExhaustive(const SeqItem& read, OverlapBlockList* pOBOut) const
{
	// Collect the complete set of overlaps in pOBOut
	static const AlignFlags sufPreAF(false, false, false);
	static const AlignFlags prePreAF(false, true, true);
	static const AlignFlags sufSufAF(true, false, true);
	static const AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	/*
	// Match the suffix of seq to prefixes
	hitTemplate.setRev(false, false);
	cost += OverlapAlgorithm::alignSuffixInexact(seq, pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	hitTemplate.setRev(true, false);
	cost += OverlapAlgorithm::alignSuffixInexact(complement(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	// Match the prefix of seq to suffixes
	hitTemplate.setRev(false, true);
	cost += OverlapAlgorithm::alignSuffixInexact(reverseComplement(seq), pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);

	hitTemplate.setRev(true, true);
	cost += OverlapAlgorithm::alignSuffixInexact(reverse(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);
	*/

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, pOBOut, pOBOut);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, pOBOut, pOBOut);

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, pOBOut, pOBOut);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, pOBOut, pOBOut);
}

// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
void OverlapAlgorithm::overlapReadIrreducible(const SeqItem& read, OverlapBlockList* pOBOut) const
{
	// The complete set of overlap blocks are collected in obTemp
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obTemp;

	static const AlignFlags sufPreAF(false, false, false);
	static const AlignFlags prePreAF(false, true, true);
	static const AlignFlags sufSufAF(true, false, true);
	static const AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	// Irreducible overlaps only
	WARN_ONCE("Irreducible-only assumptions: All reads are the same length")

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, &obTemp, pOBOut);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
	obTemp.clear();

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, &obTemp, pOBOut);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
}

// Write overlap blocks out to a file
void OverlapAlgorithm::writeOverlapBlocks(size_t readIdx, const OverlapBlockList* pList, std::ofstream& writer) const
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
                                      OverlapBlockList* pOBFinal)
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
		pOBList->clear();
		return;
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

