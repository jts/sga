#include "BWT.h"

//
// Construct the bwt of the string text and all associated data structures
//
BWT::BWT(std::string text)
{
	SAString sas("str0", 0, text);
	SAStringVector cycles;
	makeRotations(cycles, sas);
	construct(1, cycles);
}

//
BWT::BWT(SAStringVector v)
{
	SAStringVector cycles;
	for(SAStringVector::const_iterator iter = v.begin(); iter != v.end(); ++iter)
	{
		makeRotations(cycles, *iter);
	}
	construct(v.size(), cycles);
}


//
// Construct the bwt from the cycles table
//
void BWT::construct(int numStrings, SAStringVector& cycles)
{
	std::sort(cycles.begin(), cycles.end());
	m_suffixArray.resize(cycles.size());
	m_bwStr.resize(cycles.size());
	m_saStr.resize(cycles.size());
	m_occurances.resize(cycles.size());
	m_offset = numStrings;

	CharSet alphabetTmp;

	// Set up the bwt string and suffix array from the cycled strings
	for(size_t i = 0; i < cycles.size(); ++i)
	{
		SAString& curr = cycles[i];

		m_suffixArray[i] = curr.id;
		char c = curr.str[curr.str.length() - 1];
		m_bwStr[i] = c;
		m_saStr[i] = curr.str[0];
		alphabetTmp.insert(c);
	}

	// Set up the alphabet string
	for(CharSet::iterator iter = alphabetTmp.begin(); iter != alphabetTmp.end(); ++iter)
	{
		m_alphabet.append(1, *iter);
	}

	// Naive algorithms to calculate O(a,i)
	for(size_t i = 0; i < m_bwStr.size(); ++i)
	{
		for(size_t j = 0; j < m_alphabet.size(); ++j)
		{
			int count = countOcc(m_bwStr, i, m_alphabet[j]);
			m_occurances.set(m_alphabet[j], i, count);
		}
	}

	// Naive algorithm to calculate C(a)
	for(size_t j = 0; j < m_alphabet.size(); ++j)
	{
		int count = 0;
		for(size_t i = m_offset; i < m_saStr.size(); ++i)
		{
			if(m_saStr[i] < m_alphabet[j])
			{
				++count;
			}
		}
		m_predMap[m_alphabet[j]] = count;
	}

	std::cout << "SuffixArray: ";
	std::copy(m_suffixArray.begin(), m_suffixArray.end(), std::ostream_iterator<Ident>(std::cout, ","));
	std::cout << "\nBWT String:  " << m_bwStr << "\n";
	std::cout << "\nPred Map:\n";
	printMap(m_predMap);
	std::cout << "\nO(a, i): \n";
	m_occurances.print();

	std::cout << "\nTable:\n";
	for(size_t i = 0; i < m_suffixArray.size(); ++i)
	{
		//std::cout << m_suffixArray[i] << "\t" << m_saStr[i] << "\t" << m_bwStr[i] << "\n";
		std::cout << i << "\t" << cycles[i] << "\n";
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
	int r_lower = m_predMap[curr] + m_offset;
	int r_upper = r_lower + m_occurances.get(curr, m_bwStr.size() - 1) - 1;
	--j;
	std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
	for(;j >= 0; --j)
	{
		curr = w[j];
		std::cout << "RL = C(" << curr << ") + O(" << curr << ", " << r_lower - 1 << ") + " << m_offset << "\n"; 
		std::cout << "RU = C(" << curr << ") + O(" << curr << ", " << r_upper << ")\n";
		std::cout << "RL = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_lower - 1) << " + " << m_offset << "\n"; 
		std::cout << "RU = " << m_predMap[curr] << " + " << m_occurances.get(curr, r_upper) << "\n"; 

		r_lower = m_predMap[curr] + m_occurances.get(curr, r_lower - 1) + m_offset;
		r_upper = m_predMap[curr] + m_occurances.get(curr, r_upper) + m_offset - 1;

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
	int r_lower = m_predMap[curr] + m_offset;
	int r_upper = r_lower + m_occurances.get(curr, m_bwStr.size() - 1) - 1;
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

		r_lower = m_predMap[curr] + m_occurances.get(curr, r_lower - 1) + m_offset;
		r_upper = m_predMap[curr] + m_occurances.get(curr, r_upper) + m_offset - 1;

		int overlapLen = len - j;
		if(overlapLen > minOverlap)
		{
			std::cout << "Found overlap of len " << overlapLen << " to: \n";
			for(int i = r_lower; i <= r_upper; ++i)
			{
				Ident id = m_suffixArray[i];
				std::cout << "\t" << id << "\n";
			}
		}
	}
}


//
// Count the number of times c appears in string s[0, i]
//
int BWT::countOcc(std::string s, int i, char c)
{
	int count = 0;
	for(int j = 0; j <= i; ++j)
	{
		if(s[j] == c)
			++count;
	}
	return count;
}


