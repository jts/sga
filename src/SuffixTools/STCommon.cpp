//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// STCommon.cpp - Implementation of basic data structures
//
#include "STCommon.h"

//
// SAId
//

// Constructors
SAID::SAID(uint64_t i)
{
	setID(i);
}

//
SAID::SAID(uint64_t i, uint64_t p)
{
	setID(i);
	setPos(p);
}

// Set the id
void SAID::setID(uint64_t i)
{
	// Clear the HIGH bits by ANDing with the low mask
	m_val &= LOW_MASK;
	
	// Shift the new position into place and set the new value
	i <<= POS_BITS;
	m_val |= i;
}

// Set the position
void SAID::setPos(uint64_t i)
{
	// Clear the LOW bits by anding with the high mask
	m_val &= HIGH_MASK;

	// Set the new value
	m_val |= i;
}

// Get the ID
uint64_t SAID::getID() const
{
	return (m_val & HIGH_MASK) >> POS_BITS;
}

// Get the position
uint64_t SAID::getPos() const
{
	return (m_val & LOW_MASK);
}


// Output
std::ostream& operator<<(std::ostream& out, const SAID& s)
{
	out << s.getID() << "," << s.getPos();
	return out;
}



//
// GSuffix
//

// output operator for a GSuffix
std::ostream& operator<<(std::ostream& out, const GSuffix& gs)
{
	out << gs.m_label << "," << gs.m_idx;
	return out;
}

//
// SuffixString
//

// Compare 2 SuffixStrings
int operator<(const SuffixString& o1, const SuffixString& o2)
{
	return o1.str < o2.str;
}

// Output a SuffixString
std::ostream& operator<<(std::ostream& out, const SuffixString& s)
{
	out << s.id << "\t" << s.str;
	return out;
}

//
// AlphaCount
//

// Constructor
AlphaCount::AlphaCount()
{
	memset(m_counts, 0, ALPHABET_SIZE * sizeof(BaseCount));
}

// Set the value for a given base
void AlphaCount::set(char b, BaseCount v)
{
	m_counts[base2Idx(b)] = v;
}

// Increment the value for a given base
void AlphaCount::increment(char b)
{
	m_counts[base2Idx(b)]++;
};


// Get the value for a given base
BaseCount AlphaCount::get(char b) const
{
	return m_counts[base2Idx(b)];
};

std::ostream& operator<<(std::ostream& out, const AlphaCount& ac)
{
	std::copy(ac.m_counts, ac.m_counts+5, std::ostream_iterator<BaseCount>(out, " "));
	return out;
}


//
// Global functions
// 

// Convert a base to an index
AIdx base2Idx(char b)
{
	switch(b)
	{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case '$':
			return 4;
		default:
			assert(false);
	}
}

// Print a map
template<class K, class V>
void printMap(const std::map<K,V>& m)
{
	for(typename std::map<K,V>::const_iterator iter = m.begin(); iter != m.end(); ++iter)
	{
		std::cout << iter->first << "\t" << iter->second << "\n";
	}
}

// Print a vector
template<class T>
void printVector(const std::vector<T>& v)
{
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, "\n"));
}

