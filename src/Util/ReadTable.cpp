//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#include <iostream>
#include "ReadTable.h"
#include "SeqReader.h"

// Read the sequences from a file
ReadTable::ReadTable(std::string filename)
{
	SeqReader reader(filename);
	SeqItem si;
	while(reader.get(si))
	{
		std::cout << "Read sequence with ID: " << si.id << " and seq: " << si.seq << "\n";
		addRead(si);
	}
	std::cerr << "Read "<< getCount() << " sequences\n";
}

//
void ReadTable::addRead(const SeqItem& r)
{
	m_table.push_back(r);
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
