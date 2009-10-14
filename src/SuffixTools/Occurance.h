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
		Occurance() : m_sampleRate(1) {}

		// Initialize the counts from the bwt string b
		void initialize(const BWStr& bwStr, int sampleRate);

		//
		inline const AlphaCount get(const BWStr& bwStr, size_t idx) const
		{
			// Quick path
			if(idx % m_sampleRate == 0)
				return m_values[idx / m_sampleRate];

			// Calculate the nearest sample to this index
			size_t lower_idx = idx / m_sampleRate;
			size_t upper_idx = lower_idx + 1;
			size_t lower_start = lower_idx * m_sampleRate;
			size_t upper_start = upper_idx * m_sampleRate;

			AlphaCount sum;

			// Choose the closest index or force the choice to lower_idx is the upper_idx is invalid
			if((idx - lower_start < upper_start - idx) || upper_idx == m_values.size())
			{
				for(size_t j = lower_start + 1; j <= idx; ++j)
					sum.increment(bwStr[j]);
				return m_values[lower_idx] + sum;
			}
			else
			{
				for(size_t j = idx + 1; j <= upper_start; ++j)
					sum.increment(bwStr[j]);
				return m_values[upper_idx] - sum;
			}
		}
		
		inline BaseCount get(const BWStr& bwStr, char a, size_t idx) const
		{
			return get(bwStr, idx).get(a);
		}

		
		void set(char a, size_t i, BaseCount s);
		void print() const;
		size_t getByteSize() const;
		void validate(const BWStr& bwStr) const;

		friend std::ostream& operator<<(std::ostream& out, const Occurance& o);
		friend std::istream& operator>>(std::istream& in, Occurance& o);

	private:
		int m_sampleRate;
		std::vector<AlphaCount> m_values;
};


#endif
