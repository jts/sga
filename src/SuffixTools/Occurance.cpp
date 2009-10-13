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
void Occurance::increment(char a, size_t i)
{
	m_values[i].increment(a);
}

//
void Occurance::set(char a, size_t i, BaseCount s)
{
	m_values[i].set(a, s);
}

//
size_t Occurance::getByteSize() const
{
	return m_values.size() * sizeof(AlphaCount);
}

std::ostream& operator<<(std::ostream& out, const Occurance& o)
{
	out << o.m_sampleRate << "\n";
	out << o.m_values.size() << "\n";
	for(size_t i = 0; i < o.m_values.size(); ++i)
		out << o.m_values[i] << "\n";
	return out;
}

std::istream& operator>>(std::istream& in, Occurance& o)
{
	in >> o.m_sampleRate;
	size_t n;
	in >> n;
	o.m_values.resize(n);
	for(size_t i = 0; i < n; ++i)
		in >> o.m_values[i];
	return in;
}


//
void Occurance::print() const
{
	for(size_t i = 0; i < m_values.size(); i++)
	{
		std::cout << m_values[i];
	}
}
