//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Occurance.h - Data structure holding the number of times
// the letter b appears in the string S before position i
//
#ifndef OCCURANCE_H
#define OCCURANCE_H
#include "STCommon.h"

class Occurance
{
	public:
		
		// Constructors
		Occurance() {}
		Occurance(BWStr b) { initialize(b); }
		Occurance(size_t n) { m_values.resize(n); }

		// Initialize the counts from the bwt string b
		void initialize(const BWStr& b);

		// 
		int get(char a, int i) const;
		void increment(char a, int i);
		void set(char a, int i, int s);
		void print() const;
		size_t getByteSize() const;

	private:
		std::vector<AlphaCount> m_values;
};


#endif
