//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWT.h - Burrows Wheeler transform of a generalized suffix array
//
#ifndef BWT_H
#define BWT_H

#include "STCommon.h"
#include "Occurance.h"
#include "InverseSuffixArray.h"

#if 0
//
// BWT
//
class BWT
{
	public:
	
		// Constructors
		BWT(std::string t);
		BWT(const StringVector& sv);

		// Exact match
		void backwardSearch(std::string w);
		void getOverlaps(std::string w, int minOverlap);

		// L[i] -> F mapping 
		size_t LF(size_t idx) const;

		// Make all the cyclic rotations of a string
		static void makeCycles(SuffixString s, SuffixStringVector* outTable);

		// Print info about the BWT, including size
		void printInfo() const;

	private:

		void sortConstruct(int numStrings, SuffixStringVector* cycled);

		Occurance m_occurance;
		AlphaCount m_predCount;
		SuffixArray m_suffixArray;
		
		// The two representitive strings in the BWT, F is the first column of the sorted "block", L is the last
		// L is the bwt string
		BWStr m_F;
		BWStr m_L;
		size_t m_numStrings;
};
#endif
#endif
