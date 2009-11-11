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

struct ReadSeed
{
	ReadSeed(size_t i, const std::string& s) : read_idx(i), seq(s), stop_pos(s.length() - 1) { isActive = true; }

	size_t read_idx;
	std::string seq;
	size_t stop_pos;
	bool isActive;
};

typedef std::list<ReadSeed> RSList;
typedef std::vector<bool> bool_vec;

namespace AssembleExact
{
	// perform an exact mode assembly over the strings
	void assemble(const unsigned int minOverlap, const BWT* pBWT, const BWT* pRevBWT, const ReadTable* pRT, const SuffixArray* pSAI);
	
	// Mark in bv the elements of pList that are identical strings to some other element that has a higher read index
	void markInitialIdentical(bool_vec& bv, const RSList* pList, const BWT* pBWT, const SuffixArray* pSAI);

	// If pos is the position in pBWT of a full-length suffix (the entire read), return its
	// index into the suffix array index
	// Invariant: pos is a full-length suffix, which is denoted by BWT[pos] == '$'
	size_t bwtPos2ReadIdx(int64_t pos, const BWT* pBWT, const SuffixArray* pSAI);

	// Extend the seed to the left by one base, if possible
	void leftExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);

	// Extend the seed to the right by one base, if possible
	void rightExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI);
};



#endif

