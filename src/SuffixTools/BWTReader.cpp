//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReader.h - Read a BWT file from disk
//
#include "BWTReader.h"
#include "BWT.h"

//
BWTReader::BWTReader(const std::string& filename) : m_stage(IOS_NONE)
{
	m_pReader = createReader(filename);
	m_stage = IOS_HEADER;
}

//
BWTReader::~BWTReader()
{
	delete m_pReader;
}

void BWTReader::read(BWT* pBWT)
{
	size_t n;
	readHeader(pBWT->m_numStrings, n);

	pBWT->m_bwStr.resize(n);
	readBWStr(pBWT->m_bwStr);

	readPred(pBWT->m_predCount);
	readOccurrence(pBWT->m_occurrence);
}

//
void BWTReader::readHeader(size_t& num_strings, size_t& num_symbols)
{
	assert(m_stage == IOS_HEADER);
	*m_pReader >> num_strings;
	*m_pReader >> num_symbols;
	m_stage = IOS_BWSTR;	
}

//
void BWTReader::readBWStr(std::string& out_str)
{
	assert(m_stage == IOS_BWSTR);
	*m_pReader >> out_str;
	m_stage = IOS_PC;
}

//
void BWTReader::readPred(AlphaCount& out_pc)
{
	assert(m_stage == IOS_PC);
	*m_pReader >> out_pc;
	m_stage = IOS_OCC;
}

//
void BWTReader::readOccurrence(Occurrence& out_occ)
{
	assert(m_stage == IOS_OCC);
	*m_pReader >> out_occ;
	m_stage = IOS_DONE;
}
