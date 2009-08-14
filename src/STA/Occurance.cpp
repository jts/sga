//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Occurance.cpp - Data structure holding the number of times
// the letter b appears in the string S before position i
//
#include "Occurance.h"

// Initialize the counts from the bwt string b
void Occurance::initialize(const BWStr& b)
{
	m_values.resize(b.size());
	
	size_t prevIdx = 0; // the index of the previous value
	for(size_t i = 0; i < b.size(); ++i)
	{
		char currB = b[i];
		for(size_t j = 0; j < ALPHABET_SIZE; ++j)
		{
			char currA = ALPHABET[j];
			size_t v = get(currA, prevIdx);

			if(currA == currB)
				v += 1;

			set(currA, i, v);
		}
		prevIdx = i;
	}
}

// 
int Occurance::get(char a, int i) const
{
	return m_values[i].get(a);
}

//
void Occurance::increment(char a, int i)
{
	m_values[i].increment(a);
}

//
void Occurance::set(char a, int i, int s)
{
	m_values[i].set(a, s);
}

//
size_t Occurance::getByteSize() const
{
	return m_values.size() * sizeof(AlphaCount);
}

//
void Occurance::print() const
{
	for(size_t i = 0; i < m_values.size(); i++)
	{
		std::cout << m_values[i];
	}
}
