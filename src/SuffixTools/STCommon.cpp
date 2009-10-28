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
SAElem::SAElem(uint64_t i)
{
	setID(i);
}

//
SAElem::SAElem(uint64_t i, uint64_t p)
{
	setID(i);
	setPos(p);
}

// Input
std::istream& operator>>(std::istream& in, SAElem& s)
{
	uint64_t i;
	uint64_t p;
	in >> i >> p;
	s.setID(i);
	s.setPos(p);
	return in;
}

// Output
std::ostream& operator<<(std::ostream& out, const SAElem& s)
{
	if(!s.isEmpty())
		out << s.getID() << " " << s.getPos();
	else
		out << -1;
	return out;
}

//
// AlphaCount
//

//
std::ostream& operator<<(std::ostream& out, const AlphaCount& ac)
{
	std::copy(ac.m_counts, ac.m_counts+ALPHABET_SIZE, std::ostream_iterator<BaseCount>(out, " "));
	return out;
}

std::istream& operator>>(std::istream& in, AlphaCount& ac)
{
	for(size_t i = 0; i < ALPHABET_SIZE; ++i)
		in >> ac.m_counts[i];
	return in;
}


//
// Global functions
// 

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

