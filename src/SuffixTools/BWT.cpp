//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWT.cpp - Burrows Wheeler transform of a generalized suffix array
//
#include "BWT.h"
#include "Timer.h"
#include <istream>

// Parse a BWT from a file
BWT::BWT(const std::string& filename)
{
	std::ifstream in(filename.c_str());
	checkFileHandle(in, filename);	
	in >> *this;
	in.close();
}

// Construct the BWT from a suffix array
BWT::BWT(const SuffixArray* pSA, const ReadTable* pRT)
{
	Timer timer("BWT Construction");
	size_t n = pSA->getSize();
	m_numStrings = pSA->getNumStrings();
	m_bwStr.resize(n);

	// Set up the bwt string and suffix array from the cycled strings
	for(size_t i = 0; i < n; ++i)
	{
		SAElem saElem = pSA->get(i);
		const SeqItem& si = pRT->getRead(saElem.getID());

		// Get the position of the start of the suffix
		uint64_t f_pos = saElem.getPos();
		uint64_t l_pos = (f_pos == 0) ? si.seq.length() : f_pos - 1;
		m_bwStr[i] = (l_pos == si.seq.length()) ? '$' : si.seq.get(l_pos);
	}

	// initialize the occurance table
	m_occurance.initialize(m_bwStr);

	// Calculate the C(a) array
	AlphaCount tmp;
	for(size_t i = 0; i < m_bwStr.size(); ++i)
	{
		tmp.increment(m_bwStr[i]);
	}

	m_predCount.set('A', 0); // A is lexographically lowest
	m_predCount.set('C', tmp.get('A'));
	m_predCount.set('G', tmp.get('A') + tmp.get('C'));
	m_predCount.set('T', tmp.get('A') + tmp.get('C') + tmp.get('G'));
}

// Compute the last to first mapping for this BWT
size_t BWT::LF(size_t idx) const
{
	return m_bwStr[idx] != '$' ? m_predCount.get(m_bwStr[idx]) + m_occurance.get(m_bwStr[idx], idx) : 0;
}

// Perform a exact search for the string w using the backwards algorithm
void BWT::backwardSearch(std::string w) const
{
	std::cout << "Searching for " << w << "\n";
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int r_lower = m_predCount.get(curr) + m_numStrings;
	int r_upper = r_lower + m_occurance.get(curr, m_bwStr.size() - 1) - 1;
	--j;
	std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		std::cout << "RL = C(" << curr << ") + O(" << curr << ", " << r_lower - 1 << ") + " << m_numStrings << "\n"; 
		std::cout << "RU = C(" << curr << ") + O(" << curr << ", " << r_upper << ")\n";
		std::cout << "RL = " << m_predCount.get(curr) << " + " << m_occurance.get(curr, r_lower - 1) << " + " << m_numStrings << "\n"; 
		std::cout << "RU = " << m_predCount.get(curr) << " + " << m_occurance.get(curr, r_upper) << "\n"; 

		r_lower = m_predCount.get(curr) + m_occurance.get(curr, r_lower - 1) + m_numStrings;
		r_upper = m_predCount.get(curr) + m_occurance.get(curr, r_upper) + m_numStrings - 1;
		std::cout << "Curr: " << curr << " Interval now: " << r_lower << "," << r_upper << "\n";
	}

	std::cout << "Interval found: " << r_lower << "," << r_upper << "\n";
}

// Perform a search for hits to read prefixes using a backward search algorithm
void BWT::getPrefixHits(std::string w, int minOverlap, bool targetRev, bool queryRev, HitData* pHits) const
{
	// Initialize the search
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int r_lower = m_predCount.get(curr) + m_numStrings;
	int r_upper = r_lower + m_occurance.get(curr, m_bwStr.size() - 1) - 1;
	--j;

	//std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		//std::cout << "RL = C(" << curr << ") + O(" << curr << ", " << r_lower - 1 << ") + " << m_numStrings << "\n"; 
		//std::cout << "RU = C(" << curr << ") + O(" << curr << ", " << r_upper << ")\n";
		//std::cout << "RL = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_lower - 1) << " + " << m_offset << "\n"; 
		//std::cout << "RU = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_upper) << "\n"; 
		r_lower = m_predCount.get(curr) + m_occurance.get(curr, r_lower - 1) + m_numStrings;
		r_upper = m_predCount.get(curr) + m_occurance.get(curr, r_upper) + m_numStrings - 1;
		//std::cout << "Curr: " << curr << " Interval now: " << r_lower << "," << r_upper << "\n";
		
		(void)targetRev;
		(void)queryRev;
		(void)pHits;
		(void)minOverlap;
		/*
		int overlapLen = len - j;
		if(overlapLen >= minOverlap)
		{
			// Create the hit
			for(int i = r_lower; i <= r_upper; ++i)
			{
				SAElem elem = m_pSuffixArray->get(i);
				
				// Only accept hits which match the prefix of another read, indicated by the suffix element having 
				// a position of 0
				if(elem.getPos() == 0)
					pHits->addHit(Hit(elem, j, overlapLen, targetRev, queryRev));
			}
		}
		*/
	}
}

// Output operator
std::ostream& operator<<(std::ostream& out, const BWT& bwt)
{

	out << bwt.m_numStrings << "\n";
	out << bwt.m_bwStr.size() << "\n";
	out << bwt.m_bwStr << "\n";
	out << bwt.m_predCount << "\n";
	out << bwt.m_occurance;
	return out;
}


// Input operator
std::istream& operator>>(std::istream& in, BWT& bwt)
{
	in >> bwt.m_numStrings;
	size_t n;
	in >> n;
	bwt.m_bwStr.resize(n);
	in >> bwt.m_bwStr;
	in >> bwt.m_predCount;
	in >> bwt.m_occurance;
	return in;
}

// write the suffix array to a file
void BWT::write(std::string& filename)
{
	std::ofstream out(filename.c_str());
	out << *this;
	out.close();
}


// Print the BWT
void BWT::print(const ReadTable* pRT, const SuffixArray* pSA) const
{
	std::cout << "i\tL(i)\tO(-,i)\tSUFF\n";
	for(size_t i = 0; i < m_bwStr.size(); ++i)
	{
		std::cout << i << "\t" << m_bwStr[i] << "\t" << m_occurance.get(i) << pSA->getSuffix(i, pRT) << "\n";
	}
}

// Print information about the BWT
void BWT::printInfo() const
{
	size_t o_size = m_occurance.getByteSize();
	size_t p_size = sizeof(m_predCount);

	size_t bwStr_size = sizeof(m_bwStr) + m_bwStr.size();
	size_t offset_size = sizeof(m_numStrings);
	size_t total_size = o_size + p_size + bwStr_size + offset_size;
	printf("BWT Size -- occ: %zu C(a): %zu L(): %zu misc: %zu TOTAL: %zu\n",
			o_size, p_size, bwStr_size, offset_size, total_size);
	printf("N: %zu Bytes per suffix: %lf\n", m_bwStr.size(), (double)total_size / m_bwStr.size());
}
