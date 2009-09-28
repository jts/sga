//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// LCPArray - Longest Common Prefix array for a suffix array
// Entry i in the array is the length of the LCP
// between SA entry i and i + 1
//
#include "LCPArray.h"

// Naive construction algorithm
LCPArray::LCPArray(SuffixArray* pSA, ReadTable* pRT)
{
	size_t saSize = pSA->getSize();
	m_data.resize(saSize - 1);

	for(size_t i = 0; i < saSize - 1; ++i)
	{
		std::string s1 = pSA->getSuffix(i, pRT);
		std::string s2 = pSA->getSuffix(i+1, pRT);
		m_data[i] = countPrefixLength(s1, s2);
	}
}

unsigned int LCPArray::get(size_t idx) const
{
	assert(idx < m_data.size());
	return m_data[idx];
}

void LCPArray::print(SuffixArray* pSA, ReadTable* pRT) const
{
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		std::cout << i << "\t" << pSA->getSuffix(i, pRT) + "$\t"  
                       << m_data[i] << "\n";
	}
	// Print the last suffix
	std::cout << m_data.size() << "\t" << pSA->getSuffix(m_data.size(), pRT) + "$\t-\n";
}

size_t LCPArray::countPrefixLength(std::string s1, std::string s2) const
{
	size_t stop = s1.size() < s2.size() ? s1.size() : s2.size();
	size_t m = 0;
	while(s1[m] == s2[m] && m < stop)
		++m;
	return m;
}
