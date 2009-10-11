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
		SuffixArray(uint64_t i, std::string t);
		SuffixArray(const StringVector& sv);
		SuffixArray(const ReadTable& sv);
		SuffixArray(const SuffixArray& a, const SuffixArray& b);

		//
		SAElem* getHead() { return &m_data[0]; }
		SAElem get(size_t idx) const { return m_data[idx]; }
		void set(size_t idx, SAElem e) { m_data[idx] = e; }
		size_t getSize() const { return m_data.size(); }
		size_t getNumStrings() const { return m_numStrings; } 
		std::string getSuffix(size_t idx, const ReadTable* pRT) const;
		size_t getSuffixLength(const ReadTable* pRT, const SAElem elem) const;

		//
		void initialize(const ReadTable& rt);
		void initialize(size_t num_suffixes, size_t num_strings);
		void validate(const ReadTable* pRT) const;
		void sort(const ReadTable* pRT);

		// Detect all the redundant strings in the data set
		SAElemPairVec detectRedundantStrings(const ReadTable* pRT) const;

		// Detect prefix/suffix overlaps in the suffix array of a minimum length
		OverlapVector extractPrefixSuffixOverlaps(int minOverlap, const ReadTable* pRT) const;

		// Remove all the suffixes from the SA that have an id in idSet
		void removeReads(const NumericIDSet& idSet);

		// Make all the cyclic rotations of a string
		static void makeCycles(SuffixString s, SuffixStringVector* outTable);

		// Operators
		friend std::ostream& operator<<(std::ostream& out, const SuffixArray& sa);
		friend std::istream& operator>>(std::istream& in, SuffixArray& sa);

		// Print funcs
		void print() const;
		void print(const ReadTable* pRT) const;

	private:
		
		// Simple construction algorithm based on sorting the cycled rotations of the strings
		void sortConstruct(int numStrings, SuffixStringVector* cycles);

		// Data members
		SAElemVector m_data;
		size_t m_numStrings;
};

#endif 
