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
	cost += BWTAlgorithms::alignSuffixInexact(seq, pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	hitTemplate.setRev(true, false);
	cost += BWTAlgorithms::alignSuffixInexact(complement(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	// Match the prefix of seq to suffixes
	hitTemplate.setRev(false, true);
	cost += BWTAlgorithms::alignSuffixInexact(reverseComplement(seq), pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);

	hitTemplate.setRev(true, true);
	cost += BWTAlgorithms::alignSuffixInexact(reverse(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);
	*/

	// Match the suffix of seq to prefixes
	BWTAlgorithms::findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, pOBOut, pOBOut);
	BWTAlgorithms::findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, pOBOut, pOBOut);

	// Match the prefix of seq to suffixes
	BWTAlgorithms::findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, pOBOut, pOBOut);
	BWTAlgorithms::findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, pOBOut, pOBOut);
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
	BWTAlgorithms::findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, &obTemp, pOBOut);
	BWTAlgorithms::findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	BWTAlgorithms::calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
	obTemp.clear();

	// Match the prefix of seq to suffixes
	BWTAlgorithms::findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, &obTemp, pOBOut);
	BWTAlgorithms::findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	BWTAlgorithms::calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
}

