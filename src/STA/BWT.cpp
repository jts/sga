#include <sstream>
#include "BWT.h"

//
// Construct the bwt of the string text and all associated data structures
//
BWT::BWT(std::string text)
{
	SuffixString ss(0, text + "$");
	SuffixStringVector cycled;
	makeCycles(ss, &cycled);
	sortConstruct(1, &cycled);
}

//
BWT::BWT(const StringVector& sv)
{
	SuffixStringVector cycled;
	for(size_t i = 0; i < sv.size(); ++i)
	{
		SuffixString suffix(i, sv[i] + "$");
		makeCycles(suffix, &cycled);
	}
	sortConstruct(sv.size(), &cycled);
}


//
// Construct the bwt from the cycles table via simple sorting
//
void BWT::sortConstruct(int numStrings, SuffixStringVector* cycles)
{
	// Resize all the arrays
	size_t n = cycles->size();

	m_suffixArray.resize(n);
	m_bwStr.resize(n);
	m_saStr.resize(n);
	m_offset = numStrings;
	
	// Sort the array
	std::sort(cycles->begin(), cycles->end());
	// Set up the bwt string and suffix array from the cycled strings
	for(size_t i = 0; i < n; ++i)
	{
		SuffixString& curr = (*cycles)[i];

		m_suffixArray[i] = curr.id;
		char c = curr.str[curr.str.length() - 1];
		
		m_bwStr[i] = c;
		m_saStr[i] = curr.str[0];
	}

	// initialize the occurance table
	m_occurance.initialize(m_bwStr);

	// Calculate the C(a) array
	
	// Calculate the number of times character x appears in the bw string
	AlphaCount tmp;
	for(size_t i = 0; i < m_bwStr.size(); ++i)
	{
		tmp.increment(m_bwStr[i]);
	}

	m_predCount.set('A', 0); // A is lexographically lowest
	m_predCount.set('C', tmp.get('A'));
	m_predCount.set('G', tmp.get('A') + tmp.get('C'));
	m_predCount.set('T', tmp.get('A') + tmp.get('C') + tmp.get('G'));

	if(0)
	{
		std::cout << "SuffixArray: ";
		std::copy(m_suffixArray.begin(), m_suffixArray.end(), std::ostream_iterator<GSuffix>(std::cout, ","));

		std::cout << "\nBWT String:  " << m_bwStr << "\n";
		std::cout << "C(a): " << m_predCount << "\n";

		std::cout << "\nTable:\n";
		std::cout << "i\tSA\tSTR\tO(a,i)\n";
		for(size_t i = 0; i < m_suffixArray.size(); ++i)
		{
			//std::cout << m_suffixArray[i] << "\t" << m_saStr[i] << "\t" << m_bwStr[i] << "\n";
			std::cout << i << "\t" << (*cycles)[i] << "\t";
			for(size_t j = 0; j < ALPHABET_SIZE; ++j)
			{
				std::cout << m_occurance.get(ALPHABET[j], i);
			}
			std::cout << "\n";
		}
	}
}

//
// Perform a exact search for the string w using the backwards algorithm
//
void BWT::backwardSearch(std::string w)
{
	std::cout << "Searching for " << w << "\n";
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int r_lower = m_predCount.get(curr) + m_offset;
	int r_upper = r_lower + m_occurance.get(curr, m_bwStr.size() - 1) - 1;
	--j;
	std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		std::cout << "RL = C(" << curr << ") + O(" << curr << ", " << r_lower - 1 << ") + " << m_offset << "\n"; 
		std::cout << "RU = C(" << curr << ") + O(" << curr << ", " << r_upper << ")\n";
		std::cout << "RL = " << m_predCount.get(curr) << " + " << m_occurance.get(curr, r_lower - 1) << " + " << m_offset << "\n"; 
		std::cout << "RU = " << m_predCount.get(curr) << " + " << m_occurance.get(curr, r_upper) << "\n"; 

		r_lower = m_predCount.get(curr) + m_occurance.get(curr, r_lower - 1) + m_offset;
		r_upper = m_predCount.get(curr) + m_occurance.get(curr, r_upper) + m_offset - 1;
		std::cout << "Curr: " << curr << " Interval now: " << r_lower << "," << r_upper << "\n";
	}

	std::cout << "Interval found: " << r_lower << "," << r_upper << "\n";
}

//
// Perform a search for overlaps
// 
void BWT::getOverlaps(std::string w, int minOverlap)
{
	std::cout << "Searching for " << w << "\n";
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	int r_lower = m_predCount.get(curr) + m_offset;
	int r_upper = r_lower + m_occurance.get(curr, m_bwStr.size() - 1) - 1;
	--j;
	//std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		
		/*
		std::cout << "RL = C(" << curr << ") + O(" << curr << ", " << r_lower - 1 << ") + " << m_offset << "\n"; 
		std::cout << "RU = C(" << curr << ") + O(" << curr << ", " << r_upper << ")\n";
		std::cout << "RL = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_lower - 1) << " + " << m_offset << "\n"; 
		std::cout << "RU = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_upper) << "\n"; 
		*/

		r_lower = m_predCount.get(curr) + m_occurance.get(curr, r_lower - 1) + m_offset;
		r_upper = m_predCount.get(curr) + m_occurance.get(curr, r_upper) + m_offset - 1;

		int overlapLen = len - j;
		if(overlapLen > minOverlap)
		{
			std::cout << "Found overlap of len " << overlapLen << " to: \n";
			for(int i = r_lower; i <= r_upper; ++i)
			{
				GSuffix id = m_suffixArray[i];
				std::cout << "\t" << id << "\n";
			}
		}
	}
}

//
// Make all the cyclic rotations of a string, placing the result in *outTable
//
void BWT::makeCycles(SuffixString s, SuffixStringVector* outTable)
{
	int l = s.str.length();
	for(int i = 0; i < l; ++i)
	{
		SuffixString r = SuffixString(s.id.getLabel(), i, s.str.substr(i, l - i) + s.str.substr(0, i));
		outTable->push_back(r);
	}
}

//
// Print information about the BWT
//
void BWT::printInfo() const
{
	size_t o_size = m_occurance.getByteSize();
	size_t p_size = sizeof(m_predCount);
	size_t sa_size = sizeof(m_suffixArray); // Base size for the vector
	for(size_t i = 0; i < m_suffixArray.size(); ++i)
	{
		sa_size += m_suffixArray[i].getByteSize();
	}

	size_t saStr_size = sizeof(m_saStr) + m_saStr.size();
	size_t bwStr_size = sizeof(m_bwStr) + m_bwStr.size();
	size_t offset_size = sizeof(m_offset);
	size_t total_size = o_size + p_size + sa_size + saStr_size + bwStr_size + offset_size;
	printf("BWT Size -- occ: %zu C(a): %zu SA: %zu F(): %zu L(): %zu misc: %zu TOTAL: %zu\n",
			o_size, p_size, sa_size, saStr_size, bwStr_size, offset_size, total_size);
	printf("N: %zu Bytes per suffix: %lf\n", m_bwStr.size(), (double)total_size / m_bwStr.size());
}
