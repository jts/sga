//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#include "ReadTable.h"

//
void ReadTable::addRead(const Read& r)
{
	m_table.push_back(r);
}

//
const Read& ReadTable::getRead(size_t idx) const
{
	assert(idx < m_table.size());
	return m_table[idx];
}

//
size_t ReadTable::getCount() const
{
	return m_table.size();
}

// 
size_t ReadTable::getSumLengths() const
{
	size_t sum = 0;
	for(size_t i = 0; i < m_table.size(); ++i)
	{
		sum += m_table[i].seq.size();
	}
	return sum;
}
