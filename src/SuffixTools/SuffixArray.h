//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixArray - Generalized suffix array
//
#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H
#include "STCommon.h"
#include "ReadTable.h"

class LCPArray;

class SuffixArray
{
	public:
		
		//
		SuffixArray() {}
		SuffixArray(const ReadTable* pRT);

		// Construction/Validation functions
		void initialize(const ReadTable& rt);
		void initialize(size_t num_suffixes, size_t num_strings);
		void validate(const ReadTable* pRT) const;
		void sort(const ReadTable* pRT);

		// Detect prefix/suffix overlaps in the suffix array of a minimum length
		OverlapVector extractPrefixSuffixOverlaps(int minOverlap, const ReadTable* pRT) const;

		// Detect all the redundant strings in the data set
		SAElemPairVec detectRedundantStrings(const ReadTable* pRT) const;

		// Remove all the suffixes from the SA that have an id in idSet
		void removeReads(const NumericIDSet& idSet);

		// Simple accessors
		inline const SAElem& get(size_t idx) const { return m_data[idx]; }
		inline void set(size_t idx, SAElem e) { m_data[idx] = e; }
		size_t getSize() const { return m_data.size(); }
		size_t getNumStrings() const { return m_numStrings; } 
		std::string getSuffix(size_t idx, const ReadTable* pRT) const;
		size_t getSuffixLength(const ReadTable* pRT, const SAElem elem) const;
		

		// Operators
		friend std::ostream& operator<<(std::ostream& out, const SuffixArray& sa);
		friend std::istream& operator>>(std::istream& in, SuffixArray& sa);

		// Print funcs
		void print() const;
		void print(const ReadTable* pRT) const;

		// friends
		friend void saca_induced_copying(SuffixArray* pSA, const ReadTable* pRT);


	private:
		
		// Data members
		SAElemVector m_data;
		size_t m_numStrings;
};

#endif 
