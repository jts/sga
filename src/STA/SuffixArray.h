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

class SuffixArray
{
	public:
		
		//
		SuffixArray(uint64_t i, std::string t);
		SuffixArray(const StringVector& sv);
		SuffixArray(const ReadTable& sv);
		SuffixArray(const SuffixArray& a, const SuffixArray& b);

		//
		SAID get(size_t idx) const { return m_data[idx]; }
		char getF(size_t idx) const { return m_F[idx]; }
		size_t getSize() const { return m_data.size(); }
		size_t getNumStrings() const { return m_numStrings; } 

		// Validate the suffix array
		void validate(const ReadTable* pRT) const;

		// Make all the cyclic rotations of a string
		static void makeCycles(SuffixString s, SuffixStringVector* outTable);

		// Print funcs
		void print() const;
		void print(const ReadTable* pRT) const;

	private:
		
		void sortConstruct(int numStrings, SuffixStringVector* cycles);

		std::vector<SAID> m_data;
		BWStr m_F; // m_F[i] is the first character in the suffix at i
		size_t m_numStrings;
};

#endif 
