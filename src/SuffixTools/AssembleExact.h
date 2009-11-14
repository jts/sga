//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// AssembleExact - Assembly algorithm for exact sequences using a BWT
//
#ifndef ASSEMBLE_EXACT_H
#define ASSEMBLE_EXACT_H
#include "BWT.h"
#include "ReadTable.h"
#include "BWTAlgorithms.h"
#include <vector>
#include <list>

struct OverlapBlock
{
	OverlapBlock(BWTIntervalPair r, int ol) : ranges(r), overlapLen(ol) {}
	BWTIntervalPair ranges;
	int overlapLen;
};

struct ReadExt
{
	ReadExt(size_t i, const std::string& s) : read_idx(i), seq(s), w(s), isActive(true), isContained(false) {}

	size_t read_idx;
	std::string seq;
	std::string w;
	std::string preseq;
	std::string postseq;
	bool isActive;
	bool isContained;
};

typedef std::list<size_t> IndexList;
typedef std::list<BWTIntervalPair> IntervalPairList;
typedef std::list<OverlapBlock> OverlapBlockList;
typedef std::vector<ReadExt> REVec;
typedef std::vector<bool> bool_vec;

namespace AssembleExact
{
	// perform an exact mode assembly over the strings
	void assemble(const unsigned int minOverlap, const BWT* pBWT, const BWT* pRevBWT, const ReadTable* pRT, const SuffixArray* pSAI);
	
	// Mark in bv the elements of pList that are identical strings to some other element that has a higher read index
	void markInitialIdentical(REVec* pREVec, const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);

	// If pos is the position in pBWT of a full-length suffix (the entire read), return its
	// index into the suffix array index
	// Invariant: pos is a full-length suffix, which is denoted by BWT[pos] == '$'
	size_t bwtPos2ReadIdx(int64_t pos, const BWT* pBWT, const SuffixArray* pSAI);

	/*
	// Extend the seed to the left by one base, if possible
	void leftExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);

	// Extend the seed to the right by one base, if possible
	void rightExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);
	*/

	// Extend the sequence by finding a consensus right-ward extension string
	void rightExtendSG(const unsigned int minOverlap, ReadExt& elem, REVec* pREVec, 
                       const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);


	// Mark in bv the read index for the head of each interval in pList as being contained 
	void markHeadContained(bool_vec& bv, size_t int_idx, size_t read_idx, const IntervalPairList* pList, const BWT* pBWT, const SuffixArray* pSAI);


	// Collect the number of extensions for A,C,G,T of such proper prefix that matches a suffix
	// of w of at least length minOverlap. If any complete containments are found, the interval 
	// of the contained strings is placed in pIList
	AlphaCount collectExtensions(const unsigned int minOverlap, const std::string& w, 
                                 const BWT* pBWT, const BWT* pRevBWT, IntervalPairList* pIList);

	// Find the minimal extension of w described by the overlapping reads of w
	// or return the empty string if none exist.
	std::string findExtension(const unsigned int minOverlap, const std::string& w, 
	                          const BWT* pBWT, const BWT* pRevBWT, BWTIntervalPair& irr_range);


	// Get the irreducible overlaps for a given read
	void findIrreducibleOverlaps(const std::string& w, const BWT* pBWT, const BWT* pRevBWT,
	                             int minOverlap, Hit& hitTemplate, HitVector* pHits);

	// Extend each block in obl until all the irreducible overlaps have been found. 
	void processIrreducibleBlocks(OverlapBlockList& obl, const size_t qlen, const BWT* pRevBWT, Hit& hitTemplate, HitVector* pHits);

	// Return the counts of the bases between the lower and upper interval in pBWT
	AlphaCount getExtCount(BWTInterval& interval, const BWT* pBWT);

	// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
	void updateOverlapBlockRangesRight(OverlapBlockList& obList, char b, const BWT* pRevBWT);
};



#endif
