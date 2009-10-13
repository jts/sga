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
		Occurance() : m_sampleRate(1) {}
		Occurance(BWStr b) : m_sampleRate(1) { initialize(b); }
		Occurance(size_t n) : m_sampleRate(1) { m_values.resize(n); }

		// Initialize the counts from the bwt string b
		void initialize(const BWStr& b);

		// 
		inline AlphaCount get(size_t i) const
		{
			return m_values[i];
		}
		
		inline BaseCount get(char a, size_t i) const
		{
			return m_values[i].get(a);
		}

		void increment(char a, size_t i);
		void set(char a, size_t i, BaseCount s);
		void print() const;
		size_t getByteSize() const;

		friend std::ostream& operator<<(std::ostream& out, const Occurance& o);
		friend std::istream& operator>>(std::istream& in, Occurance& o);

	private:
		size_t m_sampleRate;
		std::vector<AlphaCount> m_values;
};


#endif
