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
void RLBWT::initializeFMIndex(int sampleRate)
{
    WARN_ONCE("Marker offset and position can be computed from AlphaCount?");

    // initialize the marker vector
    int num_samples = (m_numSymbols % sampleRate == 0) ? (m_numSymbols / sampleRate) : (m_numSymbols / sampleRate + 1);
    m_markers.resize(num_samples);

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols
    // however since a run may not end exactly on sampleRate boundaries,
    // we place a marker just after the run crossing the sampleRate boundary has 
    // ended. 
    
#if 0
    size_t curr_marker_index = 0;
    size_t next_marker = sampleRate;
    size_t running_total = 0;
    for(size_t i = 0; i < m_rlString.size(); ++i)
    {
        while(running_total >= next_marker)
        {
            // Place markers
            size_t expected_marker_pos = (curr_marker_index + 1) * sampleRate;
            assert(expected_marker_pos >= running_total);
            RLMarker marker;
            marker.unitIndex = i;
            marker.offset

            next_marker += sampleRate;
        }

        // Advance the running total
    }
#endif
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
    size_t m_size = m_markers.capacity() * sizeof(RLMarker);
    size_t bwStr_size = m_rlString.capacity() * sizeof(RLUnit);
    size_t other_size = sizeof(*this);
    size_t total_size = m_size + bwStr_size + other_size;

    double total_mb = ((double)total_size / (double)(1024 * 1024));
    
    printf("RLBWT contains %zu symbols in %zu runs (%1.4lf symbols per run)\n", m_numSymbols, m_rlString.size(), (double)m_numSymbols / m_rlString.size());
    printf("RLBWT Size -- Markers: %zu Str: %zu Misc: %zu TOTAL: %zu (%lf MB)\n",
            m_size, bwStr_size, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n", m_numSymbols, (double)total_size / m_numSymbols);
    
}
