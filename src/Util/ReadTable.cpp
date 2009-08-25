//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#include <iostream>
#include <algorithm>
#include "ReadTable.h"
#include "SeqReader.h"

// Read the sequences from a file
ReadTable::ReadTable(std::string filename)
{
	SeqReader reader(filename);
	SeqItem si;
	while(reader.get(si))
	{
		addRead(si);
	}
	std::cerr << "Read "<< getCount() << " sequences\n";
}

// Populate this read table with the reverse reads from pRT
void ReadTable::initializeReverse(const ReadTable* pRT)
{
	size_t numReads = pRT->getCount();
	m_table.reserve(numReads);
	for(size_t i = 0; i < numReads; ++i)
	{
		SeqItem read = pRT->getRead(i);
		read.seq = reverse(read.seq);
		addRead(read);
	}
}

//
void ReadTable::addRead(const SeqItem& r)
{
	m_table.push_back(r);
}

//
size_t ReadTable::getReadLength(size_t idx) const
{
	return m_table[idx].seq.length();
}

//
const SeqItem& ReadTable::getRead(size_t idx) const
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

//
std::ostream& operator<<(std::ostream& out, const ReadTable& rt)
{
	size_t numReads = rt.getCount();
	for(size_t i = 0; i < numReads; ++i)
	{
		const SeqItem& read = rt.getRead(i);
		out << read.id << "\t" << read.seq << "\n";
	}
	return out;
}
