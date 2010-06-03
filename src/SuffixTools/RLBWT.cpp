//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RLBWT - Run-length encoded Burrows Wheeler transform
//
#include "RLBWT.h"
#include "Timer.h"
#include "BWTReader.h"
#include "BWTWriter.h"
#include <istream>
#include <queue>
#include <inttypes.h>

// macros
#define OCC(c,i) m_occurrence.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

// Parse a BWT from a file
RLBWT::RLBWT(const std::string& filename) : m_numStrings(0), m_numSymbols(0)
{
    BWTReader reader(filename);
    reader.read(this);
}

void RLBWT::append(char b)
{
    bool increment = false;
    if(!m_rlString.empty())
    {
        RLUnit& lastUnit = m_rlString.back();
        if(lastUnit.getChar() == b && !lastUnit.isFull())
        {
            lastUnit.incrementCount();
            increment = true;
        }
    }

    if(!increment)
    {
        // Must add a new unit to the string
        RLUnit unit(b);
        m_rlString.push_back(unit);
    }
    ++m_numSymbols;
}

// Fill in the FM-index data structures
void RLBWT::initializeFMIndex()
{
    // initialize the occurance table
    m_occurrence.initialize(m_bwStr, DEFAULT_SAMPLE_RATE);

    // Calculate the C(a) array
    
    // Calculate the total number of occurances of each character in the BW str
    AlphaCount tmp;
    for(size_t i = 0; i < m_bwStr.length(); ++i)
    {
        tmp.increment(m_bwStr.get(i));
    }

    m_predCount.set('$', 0);
    m_predCount.set('A', tmp.get('$')); 
    m_predCount.set('C', m_predCount.get('A') + tmp.get('A'));
    m_predCount.set('G', m_predCount.get('C') + tmp.get('C'));
    m_predCount.set('T', m_predCount.get('G') + tmp.get('G'));
}

//
void RLBWT::validate() const
{
    std::cerr << "Warning BWT validation is turned on\n";
    m_occurrence.validate(m_bwStr);
}

// write the BWT to a file
void RLBWT::write(const std::string& filename)
{
    BWTWriter writer(filename);
    writer.write(this);
}

// Print the BWT
void RLBWT::print(const ReadTable* pRT, const SuffixArray* pSA) const
{
    std::cout << "i\tL(i)\tF(i)\tO(-,i)\tSUFF\n";
    for(size_t i = 0; i < m_bwStr.length(); ++i)
    {
        assert(getF(i) == pSA->getSuffix(i, pRT)[0]);
        std::cout << i << "\t" << m_bwStr.get(i) << "\t" << getF(i) << "\t" << m_occurrence.get(m_bwStr, i) << pSA->getSuffix(i, pRT) << "\n";
    }
}

// Print information about the BWT
void RLBWT::printInfo() const
{
    size_t o_size = m_occurrence.getByteSize();
    size_t p_size = sizeof(m_predCount);

    size_t bwStr_size = m_bwStr.getMemSize();
    size_t offset_size = sizeof(m_numStrings);
    size_t total_size = o_size + p_size + bwStr_size + offset_size;
    double total_mb = ((double)total_size / (double)(1024 * 1024));
    printf("BWT Size -- OCC: %zu C: %zu Str: %zu Misc: %zu TOTAL: %zu (%lf MB)\n",
            o_size, p_size, bwStr_size, offset_size, total_size, total_mb);
    printf("N: %zu Bytes per suffix: %lf\n", m_bwStr.length(), (double)total_size / m_bwStr.length());
}
