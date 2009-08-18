//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixArray - Generalized suffix array
//
#include "SuffixArray.h"
#include "InverseSuffixArray.h"

// Sort two suffixes
bool SuffixCompare::operator()(SAID x, SAID y) 
{ 
	const std::string& sx = m_pRT->getRead(x.getID()).seq;
	const std::string& sy = m_pRT->getRead(y.getID()).seq;

	bool xterminal = x.getPos() == sx.size();
	bool yterminal = y.getPos() == sy.size();
	if(xterminal && yterminal)
	{
		return x.getID() < y.getID();
	}
	else if(xterminal)
	{
		return true;
	}
	else if(yterminal)
	{
		return false;
	}
	else
	{
		std::string sfx = sx.substr(x.getPos()) + "$";
		std::string sfy = sy.substr(y.getPos()) + "$";
		if(sfx == sfy)
		{
			return x.getID() < y.getID();
		}
		else
		{
			return sfx < sfy;
		}
	}
}


//
// Construct the suffix array for the string
// 
SuffixArray::SuffixArray(uint64_t i, std::string text)
{
	SuffixString ss(i, text + "$");
	SuffixStringVector cycled;
	makeCycles(ss, &cycled);
	sortConstruct(1, &cycled);
}

//
// Construct the suffix array for a vector of strings using the trivial sort method
//
SuffixArray::SuffixArray(const StringVector& sv)
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
// Construct the suffix array for a table of strings using the trivial sort method
//

SuffixArray::SuffixArray(const ReadTable& rt)
{
	SuffixStringVector cycled;
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		const Read& r = rt.getRead(i);
		SuffixString suffix(i, r.seq + "$");
		makeCycles(suffix, &cycled);
	}
	sortConstruct(rt.getCount(), &cycled);
}

//
// Merge a new suffix array which is the merger of A and B
//
SuffixArray::SuffixArray(const SuffixArray& a, const SuffixArray& b)
{
	size_t n_a = a.getSize();
	size_t n_b = b.getSize();
	size_t ns_a = a.getNumStrings();
	size_t ns_b = b.getNumStrings();

	size_t n = n_a + n_b;
	m_numStrings = ns_a + ns_b;
	m_data.resize(n);
	m_F.resize(n);

	// Get the ISA for A and B
	InverseSuffixArray isa_a(a);
	InverseSuffixArray isa_b(b);

	size_t i = 0; // the current ith lowest (lexographically) suffix in A
	size_t j = 0; // likewise in B

	for(size_t k = 0; k < n; ++k)
	{
		//printf("Top of main loop (%zu, %zu, %zu)\n", i, j, k);

		// Compare A[i] to B[j]
		// Set up indices into the F array of A and B
		// These indices record the current characters to compare to determine the precedence
		// of the suffix A[i] and B[j]

		size_t fidx_a = i;
		size_t fidx_b = j;

		bool a_lowest = false;
		
		while(1)
		{
			// Early exit if A or B have all their suffixes inserted
			if(i == n_a)
			{
				a_lowest = false;
				break;
			}
			else if(j == n_b)
			{
				a_lowest = true;
				break;
			}

			char nextA = a.getF(fidx_a);
			char nextB = b.getF(fidx_b);

			//printf("	Top of while loop (%zu, %zu, %c, %c)\n", fidx_a, fidx_b, nextA, nextB);

			if(nextA == '$' && nextB == '$') // If both suffixes are terminal, return the one with the lexographically lower id
			{
				// Both suffixes are terminal, return the suffix
				// with the lowest id
				a_lowest = (a.get(i).getID() < b.get(j).getID()) ? true : false;
				break;
			}
			else if(nextA < nextB) // this catches the case where A == $ since $ is lower than the other chars
			{
				a_lowest = true;
				break;
			}
			else if(nextB < nextA) // likewise as above
			{
				a_lowest = false;
				break;
			}
			else
			{
				// The suffixes are non-terminal and equal at this position,
				// redo the compare at the next character in the suffix
				// Update fidx_a and fidx_b to point to the next character in the 
				// respective suffixes using the ISA
				// This is guarenteed to never wrap around since the match will terminate
				// at the last character
				SAID said_a = a.get(fidx_a); 
				SAID said_b = b.get(fidx_b); 

				fidx_a = isa_a.getRank(said_a.getID(), said_a.getPos() + 1);
				fidx_b = isa_b.getRank(said_b.getID(), said_b.getPos() + 1);
				//printf("	updated fidx to be (%zu, %zu)\n", fidx_a, fidx_b);
			}
		}

		// Add the next suffix
		if(a_lowest)
		{
			m_data[k] = a.get(i);
			m_F[k] = a.getF(i);
			++i;
		}
		else
		{
			m_data[k] = b.get(j);
			m_F[k] = b.getF(j);
			++j;
		}
	}
}

//
// Initialize a suffix array for the strings in RT
//
void SuffixArray::initialize(const ReadTable& rt)
{
	size_t n = rt.getSumLengths() + rt.getCount(); // We need room for 1 suffix per base pair + a '$' char per read
	m_data.resize(n);
	m_F.resize(n);
	std::cerr << "Allocating space for " << n << " suffixes\n";

	// Fill the data table with the linear ordering of the suffixes
	size_t count = 0;
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		for(size_t j = 0; j < rt.getRead(i).seq.size() + 1; ++j)
		{
			m_data[count++] = SAID(i, j);
		}
	}
	std::cerr << "Created " << count << " suffixes\n";
}

//
// Sort an initialized suffix array in-place
// 
void SuffixArray::sort(const ReadTable* pRT)
{
	SuffixCompare compare(pRT);
	std::sort(m_data.begin(), m_data.end(), compare);
}


//
// Construct the bwt from the cycles table via simple sorting
//
void SuffixArray::sortConstruct(int numStrings, SuffixStringVector* cycles)
{
	// Resize all the arrays
	size_t n = cycles->size();

	m_data.resize(n);
	m_F.resize(n);
	m_numStrings = numStrings;
	
	// Sort the array
	std::sort(cycles->begin(), cycles->end());

	// Set up the bwt string and suffix array from the cycled strings
	for(size_t i = 0; i < n; ++i)
	{
		SuffixString& curr = (*cycles)[i];
		m_data[i] = curr.id;
		m_F[i] = curr.str[0];
	}

	if(0)
	{
		std::cout << "\nTable:\n";
		std::cout << "i\tSA\tSTR\tISA\n";
		for(size_t i = 0; i < m_data.size(); ++i)
		{
			//std::cout << m_suffixArray[i] << "\t" << m_F[i] << "\t" << m_L[i] << "\n";
			//SuffixString& ss = (*cycles)[i];
			std::cout << i << "\t" << (*cycles)[i] << "\n";
		}
	}
}

//
// Validate the suffix array using the read table
//
void SuffixArray::validate(const ReadTable* pRT) const
{
	size_t maxIdx = pRT->getCount();
	size_t n = m_data.size();

	// Compute the ISA
	InverseSuffixArray isa(*this);

	// Validate the ISA is a permutation of 1..n, this implies that the id,pos pairs of the SA are valid
	isa.validate();

	// Ensure that the suffix at pos i is lexographically lower than the suffix at i + 1 using the full string
	for(size_t i = 0; i < n - 1; ++i)
	{
		SAID id1 = m_data[i];
		SAID id2 = m_data[i+1];
		assert(id1.getID() < maxIdx);
		std::string suffix1 = pRT->getRead(id1.getID()).seq.substr(id1.getPos()) + "$";
		std::string suffix2 = pRT->getRead(id2.getID()).seq.substr(id2.getPos()) + "$";
		bool suffixValidated = true;

		if(suffix1 == suffix2)
		{
			suffixValidated = id1.getID() < id2.getID();
		}
		else
		{
			suffixValidated = suffix1 < suffix2;
		}

		if(!suffixValidated)
		{
			std::cerr << "Validation failure: " << suffix1 << " is not less than " << suffix2
						<< " ids: " << id1.getID() << "," << id2.getID() << "\n";
			assert(suffix1 < suffix2);
		}
	}
}

//
// Print the suffix array
//
void SuffixArray::print(const ReadTable* pRT) const
{
	std::cout << "i\tSA(i)\tF(i)\n";
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		SAID id1 = m_data[i];
		std::string suffix1 = pRT->getRead(id1.getID()).seq.substr(id1.getPos()) + "$";
		std::cout << i << "\t" << id1 << "\t" << m_F[i] << "\t" << suffix1 << "\n";
	}
}

//
// Print the suffix array
//
void SuffixArray::print() const
{
	std::cout << "i\tSA(i)\tF(i)\n";
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		std::cout << i << "\t" << m_data[i] << "\t" << m_F[i] << "\n";
	}
}

//
// Make all the cyclic rotations of a string, placing the result in *outTable
//
void SuffixArray::makeCycles(SuffixString s, SuffixStringVector* outTable)
{
	int l = s.str.length();
	for(int i = 0; i < l; ++i)
	{
		SuffixString r = SuffixString(s.id.getID(), i, s.str.substr(i, l - i) + s.str.substr(0, i));
		outTable->push_back(r);
	}
}


