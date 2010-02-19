//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Pileup - An array of strings containing all the seen
// bases for a given region/read
//
#include "Pileup.h"
#include "Alphabet.h"
#include <math.h>

//
void Pileup::add(char b)
{
	m_data.push_back(PUElem(b, log(DEFAULT_PROB)));
}

//
void Pileup::add(char b, double lp)
{
	m_data.push_back(PUElem(b, lp));
}

// Calculate the consensus base at this position
// using a simple model where all bases are treated equally
char Pileup::calculateSimpleConsensus() const
{
	assert(!m_data.empty());
	AlphaCount ac;
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		ac.increment(m_data[i].base);
	}
	return ac.getMaxBase();
}

//
char Pileup::getCount(char base) const
{
	AlphaCount ac;
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		ac.increment(m_data[i].base);
	}
	return ac.get(base);
}

//
char Pileup::getBase(size_t idx) const
{
	assert(idx < m_data.size());
	return m_data[idx].base;
}

//
size_t Pileup::getDepth() const
{
	return m_data.size();
}

//
std::string Pileup::toStr() const
{
	std::string out;
	out.reserve(m_data.size());
	for(size_t i = 0; i < m_data.size(); ++i)
		out.append(1, m_data[i].base);
	return out;
}
