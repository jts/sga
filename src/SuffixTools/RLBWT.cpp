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
RLBWT::RLBWT(const std::string& filename) : m_numStrings(0), m_numSymbols(0), m_sampleRate(0)
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
    m_sampleRate = sampleRate;
    m_shiftValue = Occurrence::calculateShiftValue(m_sampleRate);

    // initialize the marker vector
    size_t num_samples = (m_numSymbols % m_sampleRate == 0) ? (m_numSymbols / m_sampleRate) : (m_numSymbols / m_sampleRate + 1);
    m_markers.resize(num_samples);

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols however since a run may
    // not end exactly on sampleRate boundaries, we place the markers AFTER
    // the run crossing the boundary ends
    size_t curr_marker_index = 0;

    size_t next_marker = m_sampleRate;
    size_t running_total = 0;
    AlphaCount ac;
    for(size_t i = 0; i < m_rlString.size(); ++i)
    {
        // Update the count and advance the running total
        RLUnit& unit = m_rlString[i];
        char symbol = unit.getChar();
        uint8_t run_len = unit.getCount();
        ac.add(symbol, run_len);

        running_total += run_len;

        // If this is the last symbol, place a final marker at the end of the data
        bool place_last_marker = (i == m_rlString.size() - 1) && curr_marker_index < num_samples;
        while(running_total >= next_marker || place_last_marker)
        {
            // Place markers
            size_t expected_marker_pos = (curr_marker_index + 1) * m_sampleRate;

            //int diff = running_total - expected_marker_pos;
            //printf("Placing marker at index %zu, expected pos %zu, actual pos %zu, diff %d\n", 
            //        curr_marker_index, expected_marker_pos, running_total, diff);

            // Sanity checks
            // The marker position should always be less than the running total unless 
            // the number of symbols is smaller than the sample rate
            assert(expected_marker_pos <= running_total || place_last_marker);
            assert((running_total - expected_marker_pos) <= FULL_COUNT || place_last_marker);
            assert(curr_marker_index < num_samples);
            assert(ac.getSum() == running_total);

            RLMarker& marker = m_markers[curr_marker_index];
            marker.unitIndex = i + 1;
            marker.counts = ac;
            next_marker += m_sampleRate;
            ++curr_marker_index;
            place_last_marker = false;
        }        
    }

    assert(curr_marker_index = num_samples);
    m_predCount = ac;
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
    printf("Debug -- size of marker: %zu size of alphacount: %zu\n", sizeof(RLMarker), sizeof(AlphaCount));    
}
