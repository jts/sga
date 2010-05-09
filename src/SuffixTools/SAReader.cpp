//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SAReader.h - Read a suffix array file from disk
//
#include "SAReader.h"
#include "SuffixArray.h"

//
SAReader::SAReader(const std::string& filename) : m_stage(SAIOS_NONE)
{
	m_pReader = createReader(filename);
	m_stage = SAIOS_HEADER;
}

//
SAReader::~SAReader()
{
	delete m_pReader;
}

void SAReader::read(SuffixArray* pSA)
{
	size_t num_elems;
	readHeader(pSA->m_numStrings, num_elems);
	pSA->m_data.reserve(num_elems);
	
	readElems(pSA->m_data);
}

//
void SAReader::readHeader(size_t& num_strings, size_t& num_elems)
{
	assert(m_stage == SAIOS_HEADER);
	uint16_t magic_number;

	// Ensure the file format is sane
	*m_pReader >> magic_number;
	if(magic_number != SA_FILE_MAGIC)
	{
		std::cerr << "Suffix array file is not properly formatted, aborting\n";
		exit(EXIT_FAILURE);
	}
	*m_pReader >> num_strings;
	*m_pReader >> num_elems;

	// Explicitly set the stream
	// to the start of the elements
	m_pReader->get();
	m_stage = SAIOS_ELEM;	
}

//
void SAReader::readElems(SAElemVector& elemVector)
{
	assert(m_stage == SAIOS_ELEM);
	size_t cap = elemVector.capacity();
	size_t num_read = 0;
	SAElem e;
	while(*m_pReader >> e)
	{
		elemVector.push_back(e);
		++num_read;
	}
	assert(cap >= num_read);
	m_stage = SAIOS_DONE;
}

SAElem SAReader::readElem()
{
	assert(m_stage == SAIOS_ELEM);
	SAElem e;
	if(*m_pReader >> e)
	{
		return e;
	}
	else
	{
		std::cerr << "Error: Read past end-of-file in SAReader::readElem()\n";
		exit(EXIT_FAILURE);
	}
}

