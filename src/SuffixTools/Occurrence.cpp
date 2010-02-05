//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Occurrence.cpp - Data structure holding the number of times
// the letter b appears in the string S from S[0..i] (inclusive)
//
#include "Occurrence.h"
#include "BWT.h"

// Initialize the counts from the bwt string b
void Occurrence::initialize(const BWStr& bwStr, int sampleRate)
{
	m_sampleRate = sampleRate;
	calculateShiftValue();

	size_t l = bwStr.length();
	int num_samples = (l % m_sampleRate == 0) ? (l / m_sampleRate) : (l / m_sampleRate + 1);
	m_values.resize(num_samples);
	
	AlphaCount sum;
	for(size_t i = 0; i < l; ++i)
	{
		char currB = bwStr[i];
		sum.increment(currB);
		if(i % m_sampleRate == 0)
			m_values[i / m_sampleRate] = sum;
	}
}

// 
void Occurrence::calculateShiftValue()
{
	assert(m_sampleRate > 0);
	assert(IS_POWER_OF_2(m_sampleRate));

	// m_sampleRate is a power of 2, count what bit is set
	unsigned int v = m_sampleRate;
	unsigned int c = 0; // c accumulates the total bits set in v

	while(v != 1)
	{
		v >>= 1;
		++c;
	}
	m_shift = c;
	assert(1 << m_shift == m_sampleRate);
}

//
void Occurrence::set(char a, size_t i, BaseCount s)
{
	m_values[i].set(a, s);
}

//
size_t Occurrence::getByteSize() const
{
	return m_values.size() * sizeof(AlphaCount);
}

// Validate that the sampled occurance array is correct
void Occurrence::validate(const BWStr& bwStr) const
{
	size_t l = bwStr.length();
	AlphaCount sum;
	for(size_t i = 0; i < l; ++i)
	{
		char currB = bwStr[i];
		sum.increment(currB);
		AlphaCount calculated = get(bwStr, i);
		for(int i = 0; i < ALPHABET_SIZE; ++i)
			assert(calculated.get(ALPHABET[i]) == sum.get(ALPHABET[i]));
	}
}

std::ostream& operator<<(std::ostream& out, const Occurrence& o)
{
	out << o.m_sampleRate << "\n";
	out << o.m_values.size() << "\n";
	for(size_t i = 0; i < o.m_values.size(); ++i)
		out << o.m_values[i] << "\n";
	return out;
}

std::istream& operator>>(std::istream& in, Occurrence& o)
{
	in >> o.m_sampleRate;
	size_t n;
	in >> n;
	o.m_values.resize(n);
	for(size_t i = 0; i < n; ++i)
		in >> o.m_values[i];
	o.calculateShiftValue();
	return in;
}


//
void Occurrence::print() const
{
	for(size_t i = 0; i < m_values.size(); i++)
	{
		std::cout << m_values[i];
	}
}
