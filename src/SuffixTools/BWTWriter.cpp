//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter.h - Read a BWT file from disk
//
#include "BWTWriter.h"
#include "BWT.h"

//
BWTWriter::BWTWriter(const std::string& filename) : m_stage(IOS_NONE)
{
	m_pWriter = createWriter(filename);
	m_stage = IOS_HEADER;
}

//
BWTWriter::~BWTWriter()
{
	delete m_pWriter;
}

void BWTWriter::write(const BWT* pBWT)
{
	writeHeader(pBWT->m_numStrings, pBWT->m_bwStr.size());
	writeBWStr(pBWT->m_bwStr);
	writePred(pBWT->m_predCount);
	writeOccurrence(pBWT->m_occurrence);
}

//
void BWTWriter::writeHeader(const size_t& num_strings, const size_t& num_symbols)
{
	assert(m_stage == IOS_HEADER);
	*m_pWriter << num_strings << "\n";
	*m_pWriter << num_symbols << "\n";
	m_stage = IOS_BWSTR;	
}

//
void BWTWriter::writeBWStr(const std::string& str)
{
	assert(m_stage == IOS_BWSTR);
	*m_pWriter << str << "\n";
	m_stage = IOS_PC;
}

// Write a single character of the BWStr
// If the char is '\n' we are finished
void BWTWriter::writeBWChar(char b)
{
	m_pWriter->put(b);
	if(b == '\n')
		m_stage = IOS_PC;
}

//
void BWTWriter::writePred(const AlphaCount& pc)
{
	assert(m_stage == IOS_PC);
	*m_pWriter << pc << "\n";
	m_stage = IOS_OCC;
}

//
void BWTWriter::writeOccurrence(const Occurrence& occ)
{
	assert(m_stage == IOS_OCC);
	*m_pWriter << occ; // delibrately no newline
	m_stage = IOS_DONE;
}
