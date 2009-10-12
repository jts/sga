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
#include "LCPArray.h"
#include "bucketSort.h"
#include "SuffixCompare.h"
#include "mkqs.h"

// Construct the suffix array for the string
SuffixArray::SuffixArray(uint64_t i, std::string text)
{
	SuffixString ss(i, text + "$");
	SuffixStringVector cycled;
	makeCycles(ss, &cycled);
	sortConstruct(1, &cycled);
}

// Construct the suffix array for a vector of strings using the trivial sort method
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

// Construct the suffix array for a table of strings using the trivial sort method
SuffixArray::SuffixArray(const ReadTable& rt)
{
	SuffixStringVector cycled;
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		const SeqItem& r = rt.getRead(i);
		SuffixString suffix(i, r.seq.toString() + "$");
		makeCycles(suffix, &cycled);
	}
	sortConstruct(rt.getCount(), &cycled);
}

// Initialize a suffix array for the strings in RT
void SuffixArray::initialize(const ReadTable& rt)
{
	size_t n = rt.getSumLengths() + rt.getCount(); // We need room for 1 suffix per base pair + a '$' char per read
	m_data.resize(n);
	m_numStrings = rt.getCount();
	std::cerr << "Allocating space for " << n << " suffixes\n";

	// Fill the data table with the linear ordering of the suffixes
	size_t count = 0;
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		// + 1 below is for the empty suffix (is it actually needed?)
		for(size_t j = 0; j < rt.getRead(i).seq.length() + 1; ++j)
		{
			m_data[count++] = SAElem(i, j);
		}
	}
	std::cerr << "Created " << count << " suffixes\n";
}

void SuffixArray::initialize(size_t num_suffixes, size_t num_strings)
{
	m_data.resize(num_suffixes);
	m_numStrings = num_strings;
}

// Sort an initialized suffix array in-place
void SuffixArray::sort(const ReadTable* pRT)
{
	//print(pRT);
	//std::sort(m_data.begin(), m_data.end(), compare);
	//bucketSort(m_data.begin(), m_data.end(), compare);
	SuffixCompareRadix radix_compare(pRT);
	SuffixCompareID id_compare(pRT);
	histogramSort(&m_data[0], m_data.size(), 0, radix_compare, id_compare);

	//mkqs(&m_data[0], m_data.size(), 0, compare);
	//print(pRT);
	//validate(pRT);
	//print(pRT);
	//assert(false);
}


// Construct the bwt from the cycles table via simple sorting
void SuffixArray::sortConstruct(int numStrings, SuffixStringVector* cycles)
{
	// Resize all the arrays
	size_t n = cycles->size();

	m_data.resize(n);
	m_numStrings = numStrings;
	
	// Sort the array
	std::sort(cycles->begin(), cycles->end());

	// Set up the bwt string and suffix array from the cycled strings
	for(size_t i = 0; i < n; ++i)
	{
		SuffixString& curr = (*cycles)[i];
		m_data[i] = curr.id;
	}

	if(0)
	{
		std::cout << "\nTable:\n";
		std::cout << "i\tSA\tSTR\tISA\n";
		for(size_t i = 0; i < m_data.size(); ++i)
		{
			//std::cout << m_suffixArray[i] << "\t" << m_L[i] << "\n";
			//SuffixString& ss = (*cycles)[i];
			std::cout << i << "\t" << (*cycles)[i] << "\n";
		}
	}
}

// Detect identical reads that can be removed from the collection
SAElemPairVec SuffixArray::detectRedundantStrings(const ReadTable* pRT) const
{
	// Build the LCP array
	LCPArray* pLCP = new LCPArray(this, pRT);
	
	SAElemPairVec spv;
	size_t block_root_idx = 0; // tracks the last full-length suffix seen
	size_t idx = 1;

	while(idx < m_data.size())
	{
		// Check if the current entry is identical to the previous entry 
		// (which is identical to the entry at block_root_idx)
		const SAElem& currSAElem = get(idx);
		if(currSAElem.isFull())
		{
			if(pLCP->get(idx - 1) == pRT->getReadLength(currSAElem.getID()))
			{
				// This read is identical to the block root
				spv.push_back(std::make_pair(currSAElem, get(block_root_idx)));
			}
			else
			{
				// Not identical to previous block
				block_root_idx = idx;
			}
		}
		else
		{
			block_root_idx = idx;
		}
		++idx;
	}

	delete pLCP;
	return spv;
}

// Extract all the exact prefix/suffix matches from the suffix array.
// The algorithm iterates through the SA until it finds a full-length suffix (which is also a prefix of the string/read). 
// It then iterates backwards from that position collecting all the proper suffixes that match some portion of the prefix
OverlapVector SuffixArray::extractPrefixSuffixOverlaps(int minOverlap, const ReadTable* pRT) const
{
	// Build the LCP array
	LCPArray* pLCP = new LCPArray(this, pRT);
	
	OverlapVector ov;
	size_t i = 0;
	while(i < m_data.size())
	{
		const SAElem& iElem = get(i);
		if(iElem.isFull())
		{
			assert(i > 0);
			int j = i - 1;
			// track the minimum LCP seen in the block so far
			// this is the amount the suffix at j matches the prefix at i
			size_t min_lcp = pLCP->get(j);
			std::set<uint64_t> idSet;
			while(j >= 0 && static_cast<int>(min_lcp) >= minOverlap)
			{
				// Check if the length of the longest common prefix is equal to the length of the suffix
				// (this implies that there is a perfect match between the suffix at j and the prefix at j+1)
				// Also, check that the reads are distinct
				SAElem jElem = get(j);
				if(min_lcp == getSuffixLength(pRT, jElem) && iElem.getID() != jElem.getID())
				{
					// This is a proper prefix/suffix match
					const SeqItem& iRead = pRT->getRead(iElem.getID());
					const SeqItem& jRead = pRT->getRead(jElem.getID());
					
					// Check if we have seen a match to this id yet
					// The first hit to a read will always be the longest
					if(idSet.find(jElem.getID()) == idSet.end())
					{
						ov.push_back(Overlap(jRead.id, jElem.getPos(), jRead.seq.length() - 1, jRead.seq.length(), 
											 iRead.id, 0, min_lcp - 1, iRead.seq.length()));
						idSet.insert(jElem.getID());
					}
				}
				--j;
				if(pLCP->get(j) < min_lcp)
					min_lcp = pLCP->get(j);
			}
		}
		++i;
	}

	delete pLCP;
	return ov;
}

// Output operator
std::ostream& operator<<(std::ostream& out, const SuffixArray& sa)
{
	// Write the size and number of strings
	out << sa.m_data.size() << "\n";
	out << sa.m_numStrings << "\n";

	for(size_t i = 0; i < sa.m_data.size(); ++i)
	{
		out << sa.m_data[i] << "\n";
	}
	return out;
}


// Input operator
std::istream& operator>>(std::istream& in, SuffixArray& sa)
{
	// Read the size and number of strings
	size_t n;
	in >> n;
	in >> sa.m_numStrings;

	sa.m_data.resize(n);
	size_t i = 0;
	SAElem id;
	while(in >> id)
		sa.m_data[i++] = id;
	assert(i == n);
	return in;
}

// Validate the suffix array using the read table
void SuffixArray::validate(const ReadTable* pRT) const
{
	size_t maxIdx = pRT->getCount();
	(void)maxIdx;
	size_t n = m_data.size();

	// Exit if there is nothing to do
	if(n == 0)
		return;

	// Compute the ISA
	InverseSuffixArray isa(*this);

	// Validate the ISA is a permutation of 1..n, this implies that the id,pos pairs of the SA are valid
	isa.validate();

	size_t empty_count = 0;
	// Ensure that the suffix at pos i is lexographically lower than the suffix at i + 1 using the full string
	for(size_t i = 0; i < n - 1; ++i)
	{
		SAElem id1 = m_data[i];
		SAElem id2 = m_data[i+1];
		assert(id1.getID() < maxIdx);
		std::string suffix1 = pRT->getRead(id1.getID()).seq.getSuffixString(id1.getPos());
		std::string suffix2 = pRT->getRead(id2.getID()).seq.getSuffixString(id2.getPos());

		if(suffix1.length() == 1)
			++empty_count;

		bool suffixValidated = true;
		if(suffix1 == suffix2)
		{
			suffixValidated = pRT->getRead(id1.getID()).id < pRT->getRead(id2.getID()).id;
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

	assert(m_numStrings == empty_count);
}

// 
void SuffixArray::removeReads(const NumericIDSet& idSet)
{
	// Early exit if the idset is empty
	if(idSet.empty())
		return;

	SAElemVector newData;
	newData.reserve(m_data.size());

	for(size_t idx = 0; idx < m_data.size(); ++idx)
	{
		SAElem id = m_data[idx];
		if(idSet.find(id.getID()) == idSet.end())
		{
			// not on the delete list
			newData.push_back(id);
		}
	}

	m_data.swap(newData);
	m_numStrings -= idSet.size();
}


// Get the suffix cooresponding to idx using the read table
std::string SuffixArray::getSuffix(size_t idx, const ReadTable* pRT) const
{
	SAElem id = m_data[idx];
	return pRT->getRead(id.getID()).seq.getSuffixString(id.getPos());
}

// Return the length of the suffix corresponding to elem 
size_t SuffixArray::getSuffixLength(const ReadTable* pRT, const SAElem elem) const
{
	size_t readLength = pRT->getReadLength(elem.getID());
	return readLength - elem.getPos();
}


// Print the suffix array
void SuffixArray::print(const ReadTable* pRT) const
{
	std::cout << "i\tSA(i)\n";
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		SAElem id1 = m_data[i];
		std::string suffix = !id1.isEmpty() ? getSuffix(i, pRT) : "";
		std::cout << i << "\t" << id1 << "\t" << suffix << "\n";
	}
}

// Print the suffix array
void SuffixArray::print() const
{
	std::cout << "i\tSA(i)\n";
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		std::cout << i << "\t" << m_data[i] << "\n";
	}
}

// Make all the cyclic rotations of a string, placing the result in *outTable
void SuffixArray::makeCycles(SuffixString s, SuffixStringVector* outTable)
{
	int l = s.str.length();
	for(int i = 0; i < l; ++i)
	{
		SuffixString r = SuffixString(s.id.getID(), i, s.str.substr(i, l - i) + s.str.substr(0, i));
		outTable->push_back(r);
	}
}


