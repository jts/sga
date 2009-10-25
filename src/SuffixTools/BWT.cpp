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
#include <queue>
#include <inttypes.h>

struct PartialAlign
{
	PartialAlign(int i, int n, int64_t s, int64_t e) : idx(i), z(n), i_start(s), i_end(e) {}
	int idx; // the index of the current base being processed
	int z; // the number of mismatches that have occured so far
	int64_t i_start; // the start of the interval in the SA
	int64_t i_end; // the end of the interval in the SA
};

typedef std::queue<PartialAlign> PAQueue;

// macros
#define OCC(c,i) m_occurance.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

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
	m_occurance.initialize(m_bwStr, DEFAULT_SAMPLE_RATE);

	// Calculate the C(a) array
	
	// Calculate the total number of occurances of each character in the BW str
	AlphaCount tmp;
	for(size_t i = 0; i < m_bwStr.size(); ++i)
	{
		tmp.increment(m_bwStr[i]);
	}

	m_predCount.set('$', 0);
	m_predCount.set('A', tmp.get('$')); 
	m_predCount.set('C', m_predCount.get('A') + tmp.get('A'));
	m_predCount.set('G', m_predCount.get('C') + tmp.get('C'));
	m_predCount.set('T', m_predCount.get('G') + tmp.get('G'));
}

// Compute the last to first mapping for this BWT
size_t BWT::LF(size_t idx) const
{
	return m_bwStr[idx] != '$' ? PRED(m_bwStr[idx]) + OCC(m_bwStr[idx], idx) : 0;
}

// Perform a exact search for the string w using the backwards algorithm
void BWT::backwardSearch(std::string w) const
{
	std::cout << "Searching for " << w << "\n";
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int r_lower = PRED(curr);
	int r_upper = r_lower + OCC(curr, m_bwStr.size() - 1) - 1;
	--j;
	std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		printf("RL = C(%c) + O(%c,%d) + %zu\n", curr, curr, r_lower - 1, m_numStrings); 
		printf("RU = C(%c) + O(%c,%d)\n", curr, curr, r_upper); 
		printf("RL = %zu + %zu + %zu\n", (size_t)PRED(curr), (size_t)OCC(curr, r_lower - 1), m_numStrings); 
		printf("RU = %zu + %zu\n", (size_t)PRED(curr), (size_t)OCC(curr, r_upper)); 
		r_lower = PRED(curr) + OCC(curr, r_lower - 1);
		r_upper = PRED(curr) + OCC(curr, r_upper) - 1;
		printf("Curr: %c, Interval now: %d,%d\n", curr, r_lower, r_upper);
	}

	std::cout << "Interval found: " << r_lower << "," << r_upper << "\n";
}

// Perform a search for hits to read prefixes using a backward search algorithm
void BWT::getPrefixHits(size_t readIdx, std::string w, int minOverlap, bool targetRev, bool queryRev, HitVector* pHits) const
{
	// Initialize the search
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int64_t r_lower = PRED(curr);
	int64_t r_upper = r_lower + OCC(curr, m_bwStr.size() - 1) - 1;
	--j;
	//std::cout << "Searching for string: " << w << "\n";
	//printf("Starting point: %zu,%zu\n", r_lower, r_upper);
	//std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];

		//printf("RL = C(%c) + O(%c,%zu) + %zu\n", curr, curr, r_lower - 1, m_numStrings); 
		//printf("RU = C(%c) + O(%c,%zu)\n", curr, curr, r_upper); 
		//printf("RL = %zu + %zu + %zu\n", PRED(curr), OCC(curr, r_lower - 1), m_numStrings); 
		//printf("RU = %zu + %zu\n", PRED(curr), OCC(curr, r_upper));

		r_lower = PRED(curr) + OCC(curr, r_lower - 1);
		r_upper = PRED(curr) + OCC(curr, r_upper) - 1;

		//printf("Curr: %c, Interval now: %zu,%zu\n", curr, r_lower, r_upper);
		int overlapLen = len - j;
		if(overlapLen >= minOverlap)
		{
			// Output the hits where the suffix of w has matched a proper prefix 
			// (starting from the begining of the string) of some other string
			// These suffixes can be calculated using the fm-index like any other interval
			int64_t t_lower = PRED('$') + OCC('$', r_lower - 1);
			int64_t t_upper = PRED('$') + OCC('$', r_upper) - 1;
			for(int64_t sa_idx = t_lower; sa_idx <= t_upper; ++sa_idx)
				pHits->push_back(Hit(readIdx, sa_idx, j, overlapLen, targetRev, queryRev, 0));
		}
	}
}

// Perform a search for hits to read prefixes using a backward search algorithm
int BWT::getInexactPrefixHits(std::string w, const BWT* pRevBWT, int maxDiff, int minOverlap, size_t readIdx, bool targetRev, bool queryRev, HitVector* pHits) const
{
	int cost = 0;
	int len = w.size();
	int* pD = new int[len];
	calculateD(w, minOverlap, pRevBWT, queryRev == targetRev, pD);
	
	int j = len - 1;
	PAQueue hitQueue;
	
	// Create the initial partial hits
	for(int i = 0; i < 4; ++i)
	{
		char base = ALPHABET[i];
		int64_t r_lower = PRED(base);
		int64_t r_upper = r_lower + OCC(base, m_bwStr.size() - 1) - 1;
		if(ALPHABET[i] == w[j])
			hitQueue.push(PartialAlign(j, maxDiff, r_lower, r_upper));
		else
		{
			if(maxDiff > 0)
				hitQueue.push(PartialAlign(j, maxDiff - 1, r_lower, r_upper));
		}
	}

	while(!hitQueue.empty())
	{
		++cost;
		PartialAlign pa = hitQueue.front();
		hitQueue.pop();

		//printf("<curr hit> idx: %d nm: %d rl: %d ru: %d\n", pa.idx, pa.num_mismatches, (int)pa.i_start, (int)pa.i_end);
		j = pa.idx;
		int overlap_len = len - j;

		if(pa.z < pD[j])
		{
			continue;
		}

		// Output valid prefix matches for this hit
		if(overlap_len >= minOverlap)
		{
			int64_t t_lower = PRED('$') + OCC('$', pa.i_start - 1);
			int64_t t_upper = PRED('$') + OCC('$', pa.i_end) - 1;
			for(int64_t sa_idx = t_lower; sa_idx <= t_upper; ++sa_idx)
				pHits->push_back(Hit(readIdx, sa_idx, j, overlap_len, targetRev, queryRev, maxDiff - pa.z));
		}

		// Calculate the next partial alignments
		--j;
		if(j >= 0)
		{
			for(int i = 0; i < 4; ++i)
			{
				char base = ALPHABET[i];
				size_t pb = PRED(base);
				int64_t r_lower = pb + OCC(base, pa.i_start - 1);
				int64_t r_upper = pb + OCC(base, pa.i_end) - 1;
				if(r_lower <= r_upper)
				{
					if(ALPHABET[i] == w[j])
						hitQueue.push(PartialAlign(j, pa.z, r_lower, r_upper));
					else
						hitQueue.push(PartialAlign(j, pa.z - 1, r_lower, r_upper));
				}
			}
		}
	}
	delete [] pD;
	return cost;
}

void BWT::calculateD(std::string w, int minOverlap, const BWT* pRevBWT, bool contains_w, int* pD) const
{
	//std::cout << "D: " << w << " contains " << contains_w << "\n";
	int min_span = (contains_w) ? 1 : 0;
	size_t len = w.length();

	int64_t r_lower = 0; 
	int64_t r_upper = 0;
	int z = 0;

	//std::cout << "w: " << w << "\n   ";
	for(size_t i = 0; i < len; ++i)
		pD[i] = 0;

	for(size_t i = len - minOverlap; i < len; ++i)
	{
		char b = w[i];
		
		if(i == len - minOverlap)
		{
			r_lower = PRED(b);
			r_upper = r_lower + pRevBWT->getOcc(b, m_bwStr.size() - 1) - 1;
		}
		else
		{
			r_lower = PRED(b) + pRevBWT->getOcc(b, r_lower - 1);
			r_upper = PRED(b) + pRevBWT->getOcc(b, r_upper) - 1;
		}

		//printf("j: %zu Curr: %c, Interval now: %zu,%zu\n", j, b, r_lower, r_upper);
		int span = r_upper - r_lower + 1;

		if(span <= min_span)
		{
			r_lower = 1;
			r_upper = m_bwStr.size() - 1;
			z += 1;
		}
		pD[i] = z;
		//pD[i] = 0;
	}
	
	/*
	for(size_t i = 0; i < len; ++i)
	{
		std::cout << pD[i];
	}
	std::cout << "\n";
	*/
}



void BWT::validate() const
{
	std::cerr << "Warning BWT validation is turned on\n";
	m_occurance.validate(m_bwStr);
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
		std::cout << i << "\t" << m_bwStr[i] << "\t" << m_occurance.get(m_bwStr, i) << pSA->getSuffix(i, pRT) << "\n";
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
	double total_mb = ((double)total_size / (double)(1024 * 1024));
	printf("BWT Size -- OCC: %zu C: %zu Str: %zu Misc: %zu TOTAL: %zu (%lf MB)\n",
			o_size, p_size, bwStr_size, offset_size, total_size, total_mb);
	printf("N: %zu Bytes per suffix: %lf\n", m_bwStr.size(), (double)total_size / m_bwStr.size());
}
