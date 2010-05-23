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
    BWFlag flag;
    readHeader(pBWT->m_numStrings, n, flag);

    pBWT->m_bwStr.resize(n);
    readBWStr(pBWT->m_bwStr);
    
    // If the flag indicates the FM-index is stored in the file
    // read it, otherwise construct it from BWStr
    if(flag == BWF_HASFMI)
    {
        readPred(pBWT->m_predCount);
        readOccurrence(pBWT->m_occurrence);
    }
    else
    {
        pBWT->initializeFMIndex();
    }
}

//
void BWTReader::readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    uint16_t magic_number;

    // Ensure the file format is sane
    *m_pReader >> magic_number;
    if(magic_number != BWT_FILE_MAGIC)
    {
        std::cerr << "BWT file is not properly formatted, aborting\n";
        exit(EXIT_FAILURE);
    }
    *m_pReader >> num_strings;
    *m_pReader >> num_symbols;
    int temp;
    *m_pReader >> temp;
    flag = static_cast<BWFlag>(temp);
    // We must explicitly set the stream
    // to the start of the BWStr portion of
    // the file. The above extraction
    // does not discard the trailing newline 
    // so we do it here
    m_pReader->get();

    m_stage = IOS_BWSTR;    
}

//
void BWTReader::readBWStr(BWTString& out_str)
{
    assert(m_stage == IOS_BWSTR);
    char b;
    size_t idx = 0;
    while(1)
    {
        m_pReader->get(b);
        if(b != '\n')
            out_str.set(idx++, b);
        else
            break;
    }

    assert(idx == out_str.length());
    m_stage = IOS_PC;
}

// Read a single base from the BWStr
char BWTReader::readBWChar()
{
    assert(m_stage == IOS_BWSTR);
    char b;
    m_pReader->get(b);
    if(b == '\n')
        m_stage = IOS_PC;
    return b;
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
