//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SeqReader - Reads fasta or fastq sequence files
//
#include <iostream>
#include "SeqReader.h"
#include "Util.h"

SeqReader::SeqReader(std::string filename) : m_fileHandle(filename.c_str())
{
	checkFileHandle(m_fileHandle, filename);
}
	
SeqReader::~SeqReader()
{
	m_fileHandle.close();
}

// Extract an element from the file
// Return true if successful
bool SeqReader::get(SeqRecord& sr)
{
	RecordType rt = RT_UNKNOWN;
	std::string header;
	while(m_fileHandle.good())
	{
		getline(m_fileHandle, header);
		if(header.empty())
			continue;

		if(header[0] == '>')
		{
			rt = RT_FASTA;
			break;
		}
		else if(header[0] == '@')
		{
			rt = RT_FASTQ;
			break;
		}
	}

	if(rt == RT_UNKNOWN)
	{
		// No valid start found
		return false;
	}
	
	// Parse the rest of the record
	bool validRecord = false;
	std::string seq;
	std::string qual;

	if(rt == RT_FASTA)
	{
		std::string temp;
		while(m_fileHandle.good() && m_fileHandle.peek() != '>')
		{
			getline(m_fileHandle, temp);
			if(m_fileHandle.good() && temp.size() > 0)
				seq.append(temp);
		}

		// The record is valid if we extracted at least 1 bp for the sequence
		// at this point the eof may have been hit but we still want the data
		validRecord = seq.size() > 0; 
	}
	else if(rt == RT_FASTQ)
	{
		std::string temp;
		getline(m_fileHandle, seq);
		getline(m_fileHandle, temp); //discard
		getline(m_fileHandle, qual);

		// FASTQ is required to have 4 fields, we must not have hit the EOF by this point
		validRecord = seq.size() > 0 && seq.size() == qual.size() && !m_fileHandle.eof();
	}

	if(validRecord)
	{
		// Parse the id
		size_t endPos = std::min(header.find_first_of(' '), header.find_first_of('\t'));
		if(endPos != std::string::npos)
		{
			assert(endPos > 0);
			sr.id = header.substr(1, endPos - 1);
		}
		else
		{
			sr.id = header.substr(1);
		}
		sr.seq = seq;
		sr.qual = qual;
	}

	return validRecord;
}
