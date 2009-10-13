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
#include "SuffixArray.h"
#include "ReadTable.h"
#include "HitData.h"

//
// BWT
//
class BWT
{
	public:
	
		// Constructors
		BWT(const std::string& filename);
		BWT(const SuffixArray* pSA, const ReadTable* pRT);
			
		// Exact match
		void backwardSearch(std::string w) const;
		void getPrefixHits(std::string w, int minOverlap, bool targetRev, bool queryRev, HitData* pHits) const;

		// L[i] -> F mapping 
		size_t LF(size_t idx) const;

		// Print info about the BWT, including size
		void printInfo() const;
		void print(const ReadTable* pRT, const SuffixArray* pSA) const;

		// IO
		friend std::ostream& operator<<(std::ostream& out, const BWT& bwt);
		friend std::istream& operator>>(std::istream& in, BWT& bwt);

	private:
		BWT() {}
		Occurance m_occurance;
		AlphaCount m_predCount;
		
		// The two representitive strings in the BWT, F is the first column of the sorted "block", L is the last
		BWStr m_bwStr;
		size_t m_numStrings;
};
#endif
