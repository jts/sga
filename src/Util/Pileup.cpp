//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Pileup - An array of strings containing all the seen
// bases for a given region/read
//
#include "Pileup.h"

//
void Pileup::add(char b, double lp)
{
	m_data.push_back(PUElem(b, lp));
}

//
std::string Pileup::getStr() const
{
	std::string out;
	out.reserve(m_data.size());
	for(size_t i = 0; i < m_data.size(); ++i)
		out.append(1, m_data[i].base);
	return out;
}
