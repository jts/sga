//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Occurance.cpp - Data structure holding the number of times
// the letter b appears in the string S from S[0..i] (inclusive)
//
#include "Occurance.h"
#include "BWT.h"

// Initialize the counts from the bwt string b
void Occurance::initialize(const BWT* pBWT, int sampleRate)
{
	m_pBWT = pBWT;
	m_sampleRate = sampleRate;
	const BWStr* pBWString = m_pBWT->getBWStr();
	size_t l = pBWString->length();
	int num_samples = (l % m_sampleRate == 0) ? (l / m_sampleRate) : (l / m_sampleRate + 1);
	m_values.resize(num_samples);
	
	AlphaCount sum;
	for(size_t i = 0; i < l; ++i)
	{
		char currB = (*pBWString)[i];
		sum.increment(currB);
		if(i % m_sampleRate == 0)
			m_values[i / m_sampleRate] = sum;
	}

	std::cerr << "Warning occurance validation is turned on\n";
	validate();
}

// Calculate the nearest sampled index and compute the values from 
// the requested index to the nearest sample
const AlphaCount Occurance::get(size_t idx) const
{
	// Quick path
	if(idx % m_sampleRate == 0)
		return m_values[idx / m_sampleRate];

	// Calculate the nearest sample to this index
	const BWStr* pBWStr = m_pBWT->getBWStr();
	size_t lower_idx = idx / m_sampleRate;
	size_t upper_idx = lower_idx + 1;
	size_t lower_start = lower_idx * m_sampleRate;
	size_t upper_start = upper_idx * m_sampleRate;

	AlphaCount sum;

	// Choose the closest index or force the choice to lower_idx is the upper_idx is invalid
	if((idx - lower_start < upper_start - idx) || upper_idx == m_values.size())
	{
		for(size_t j = lower_start + 1; j <= idx; ++j)
			sum.increment((*pBWStr)[j]);
		return m_values[lower_idx] + sum;
	}
	else
	{
		for(size_t j = idx + 1; j <= upper_start; ++j)
			sum.increment((*pBWStr)[j]);
		return m_values[upper_idx] - sum;
	}
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

// Validate that the sampled occurance array is correct
void Occurance::validate() const
{
	const BWStr* pBWString = m_pBWT->getBWStr();
	size_t l = pBWString->length();
	AlphaCount sum;
	for(size_t i = 0; i < l; ++i)
	{
		char currB = (*pBWString)[i];
		sum.increment(currB);
		AlphaCount calculated = get(i);
		for(int i = 0; i < ALPHABET_SIZE; ++i)
			assert(calculated.get(ALPHABET[i]) == sum.get(ALPHABET[i]));
	}
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
