//
//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadInfoTable - A 0-indexed table of ID, length pairs
// Used to convert suffix array hits to overlaps
//
#include <iostream>
#include <algorithm>
#include "ReadInfoTable.h"
#include "SeqReader.h"

// Read the sequences from a file
ReadInfoTable::ReadInfoTable(std::string filename, size_t num_expected)
{
    if(num_expected > 0)
    {
        m_table.reserve(num_expected);
    }

    SeqReader reader(filename);
    SeqRecord sr;
    while(reader.get(sr))
        addReadInfo(ReadInfo(sr.id, sr.seq.length()));
}

// 
ReadInfoTable::~ReadInfoTable()
{
}

//
void ReadInfoTable::addReadInfo(const ReadInfo& r)
{
    m_table.push_back(r);
}

//
size_t ReadInfoTable::getReadLength(size_t idx) const
{
    assert(idx < m_table.size());
    return m_table[idx].length;
}

//
std::string ReadInfoTable::getReadID(size_t idx) const
{
    assert(idx < m_table.size());
    return m_table[idx].id;
}

//
const ReadInfo& ReadInfoTable::getReadInfo(size_t idx) const
{
    assert(idx < m_table.size());
    return m_table[idx];
}

//
size_t ReadInfoTable::getCount() const
{
    return m_table.size();
}

// 
size_t ReadInfoTable::countSumLengths() const
{
    size_t sum = 0;
    for(size_t i = 0; i < m_table.size(); ++i)
        sum += m_table[i].length;
    return sum;
}

//
void ReadInfoTable::clear()
{
    m_table.clear();
}

