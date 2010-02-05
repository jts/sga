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
#include "Occurrence.h"
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
		void getPrefixHits(size_t readIdx, std::string w, int minOverlap, bool targetRev, bool queryRev, HitVector* pHits) const;
		int getInexactPrefixHits(std::string w, const BWT* pRevBWT, int maxDiff, int minOverlap, size_t readIdx, bool targetRev, bool queryRev, HitVector* pHits) const;

		// L[i] -> F mapping 
		size_t LF(size_t idx) const;

		inline char getChar(size_t idx) const { return m_bwStr[idx]; }
		inline BaseCount getPC(char b) const { return m_predCount.get(b); }
		inline BaseCount getOcc(char b, size_t idx) const { return m_occurance.get(m_bwStr, b, idx); }
		inline AlphaCount getOccDiff(size_t idx0, size_t idx1) const { return m_occurance.getDiff(m_bwStr, idx0, idx1); }
		inline size_t getBWLen() const { return m_bwStr.length(); }

		// Return the first letter of the suffix starting at idx
		inline char getF(size_t idx) const
		{
			size_t ci = 0;
			while(ci < ALPHABET_SIZE && m_predCount.getByIdx(ci) <= idx)
				ci++;
			assert(ci != 0);
			return RANK_ALPHABET[ci - 1];
		}


		// 
		inline const BWStr* getBWStr() const
		{
			return &m_bwStr;
		}

		// Print info about the BWT, including size
		void printInfo() const;
		void print(const ReadTable* pRT, const SuffixArray* pSA) const;
		void validate() const;

		// IO
		friend std::ostream& operator<<(std::ostream& out, const BWT& bwt);
		friend std::istream& operator>>(std::istream& in, BWT& bwt);
		void write(std::string& filename);

		size_t getNumLoops() const { return m_totalLoops; }

	private:

		// calculate the lower bound of number of differences in w[0,i]
		// if contains_w is true, the string (or read) w is contained in the bwt
		// and should not be counted
		void calculateD(std::string w, int minOverlap, const BWT* pRevBWT, bool contains_w, int* pD) const;

		static const int DEFAULT_SAMPLE_RATE = 64;
		// Default constructor is not allowed
		BWT() {}

		// The O(a,i) array
		Occurrence m_occurance;

		// The C(a) array
		AlphaCount m_predCount;
		
		// The bw string
		BWStr m_bwStr;

		// The number of strings in the collection
		size_t m_numStrings;

		// profiling
		mutable size_t m_totalLoops; 
};
#endif
