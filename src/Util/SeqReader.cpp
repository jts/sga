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

SeqReader::SeqReader(std::string filename)
{
    m_pHandle = createReader(filename);
}
    
SeqReader::~SeqReader()
{
    delete m_pHandle;
}

// Extract an element from the file
// Return true if successful
bool SeqReader::get(SeqRecord& sr)
{
    static int warn_count = 0;
    const int MAX_WARN = 10;
    RecordType rt = RT_UNKNOWN;
    std::string header;
    while(m_pHandle->good())
    {
        getline(*m_pHandle, header);
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
        while(m_pHandle->good() && m_pHandle->peek() != '>')
        {
            getline(*m_pHandle, temp);
            if(m_pHandle->good() && temp.size() > 0)
                seq.append(temp);
        }

        // The record is valid if we extracted at least 1 bp for the sequence
        // at this point the eof may have been hit but we still want the data
        validRecord = seq.size() > 0; 
    }
    else if(rt == RT_FASTQ)
    {
        std::string temp;
        getline(*m_pHandle, seq);
        getline(*m_pHandle, temp); //discard
        getline(*m_pHandle, qual);

        // FASTQ is required to have 4 fields, we must not have hit the EOF by this point
        if(seq.size() != qual.size() && warn_count++ < MAX_WARN)
        {
            std::cerr << "Warning, FASTQ quality string is not the same length as the sequence string for read " << header << "\n";
            
        }
        validRecord = seq.size() > 0 && qual.size() > 0 && !m_pHandle->eof();
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
