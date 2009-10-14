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

class BWT;

class Occurance
{
	public:
		
		// Constructors
		Occurance() : m_pBWT(NULL), m_sampleRate(1) {}

		// Initialize the counts from the bwt string b
		void initialize(const BWT* pBWT, int sampleRate);

		//
		const AlphaCount get(size_t idx) const;
		inline BaseCount get(char a, size_t idx) const
		{
			return get(idx).get(a);
		}

		
		void set(char a, size_t i, BaseCount s);
		void setBWT(BWT* pBWT) { m_pBWT = pBWT; }
		void print() const;
		size_t getByteSize() const;
		void validate() const;

		friend std::ostream& operator<<(std::ostream& out, const Occurance& o);
		friend std::istream& operator>>(std::istream& in, Occurance& o);

	private:
		const BWT* m_pBWT;
		int m_sampleRate;
		std::vector<AlphaCount> m_values;
};


#endif
