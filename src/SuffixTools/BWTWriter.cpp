//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter - Write a BWT file to disk
//
#include "BWTWriter.h"
#include "SBWT.h"
#include "RLBWT.h"

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

void BWTWriter::write(const SBWT* pBWT)
{
    writeHeader(pBWT->m_numStrings, pBWT->m_bwStr.length(), BWF_HASFMI);
    writeBWStr(pBWT->m_bwStr);
    writePred(pBWT->m_predCount);
    writeOccurrence(pBWT->m_occurrence);
}

void BWTWriter::write(const RLBWT* pRLBWT)
{
    writeHeader(pRLBWT->m_numStrings, pRLBWT->m_numSymbols, BWF_NOFMI);
    size_t numRuns = pRLBWT->getNumRuns();
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = pRLBWT->m_rlString[i];
        char symbol = unit.getChar();
        size_t length = unit.getCount();
        for(size_t j = 0; j < length; ++j)
            writeBWChar(symbol);
    }
    // Finalize the string
    writeBWChar('\n');
}

//
void BWTWriter::write(const SuffixArray* pSA, const ReadTable* pRT)
{
    size_t num_symbols = pSA->getSize();
    size_t num_strings = pSA->getNumStrings();
    writeHeader(num_strings, num_symbols, BWF_NOFMI);

    for(size_t i = 0; i < num_symbols; ++i)
    {
        SAElem saElem = pSA->get(i);
        const SeqItem& si = pRT->getRead(saElem.getID());

        // Get the position of the start of the suffix
        uint64_t f_pos = saElem.getPos();
        uint64_t l_pos = (f_pos == 0) ? si.seq.length() : f_pos - 1;
        char b = (l_pos == si.seq.length()) ? '$' : si.seq.get(l_pos);
        writeBWChar(b);
    }
}

//
void BWTWriter::writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    *m_pWriter << BWT_FILE_MAGIC << "\n";
    *m_pWriter << num_strings << "\n";
    *m_pWriter << num_symbols << "\n";
    int temp = flag;
    *m_pWriter << temp << "\n";
    m_stage = IOS_BWSTR;    
}

//
void BWTWriter::writeBWStr(const BWTString& str)
{
    assert(m_stage == IOS_BWSTR);
    size_t n = str.length();
    for(size_t i = 0; i < n; ++i)
        *m_pWriter << str.get(i);
    *m_pWriter << "\n";
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
