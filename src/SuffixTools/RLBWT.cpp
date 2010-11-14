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
#include "BWTReader.h"
#include <istream>
#include <queue>
#include <inttypes.h>

// macros
#define OCC(c,i) m_occurrence.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

// Parse a BWT from a file
RLBWT::RLBWT(const std::string& filename, int sampleRate) : m_numStrings(0), m_numSymbols(0), m_sampleRate(sampleRate)
{
    IBWTReader* pReader = BWTReader::createReader(filename);
    pReader->read(this);
    initializeFMIndex();
    delete pReader;
}

//
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
    m_shiftValue = Occurrence::calculateShiftValue(m_sampleRate);

    // initialize the marker vector, we place a marker at the beginning (with no accumulated counts), every m_sampleRate
    // bases and one at the very end (with the total counts)
    size_t num_samples = (m_numSymbols % m_sampleRate == 0) ? (m_numSymbols / m_sampleRate) + 1 : (m_numSymbols / m_sampleRate) + 2;
    m_markers.resize(num_samples);

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols however since a run may
    // not end exactly on sampleRate boundaries, we place the markers AFTER
    // the run crossing the boundary ends

    // Place a blank first marker at the start of the data
    RLMarker& marker = m_markers[0];
    marker.unitIndex = 0;

    size_t curr_marker_index = 1;
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

        // If this is the last symbol, place the final marker(s)
        bool place_last_marker = (i == m_rlString.size() - 1) && curr_marker_index < num_samples;
        while(running_total >= next_marker || place_last_marker)
        {
            // Place markers
            size_t expected_marker_pos = curr_marker_index * m_sampleRate;

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
            place_last_marker = (i == m_rlString.size() - 1) && curr_marker_index < num_samples;
        }        
    }

    assert(curr_marker_index == num_samples);

    // Initialize C(a)
    m_predCount.set('$', 0);
    m_predCount.set('A', ac.get('$')); 
    m_predCount.set('C', m_predCount.get('A') + ac.get('A'));
    m_predCount.set('G', m_predCount.get('C') + ac.get('C'));
    m_predCount.set('T', m_predCount.get('G') + ac.get('G'));    
}

// Print the BWT
void RLBWT::print() const
{
    size_t numRuns = getNumRuns();
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = m_rlString[i];
        char symbol = unit.getChar();
        size_t length = unit.getCount();
        for(size_t j = 0; j < length; ++j)
            std::cout << symbol;
        std::cout << " : " << symbol << "," << length << "\n"; 
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
    
    printf("\nRLBWT info:\n");
    printf("Sample rate: %zu\n", m_sampleRate);
    printf("Contains %zu symbols in %zu runs (%1.4lf symbols per run)\n", m_numSymbols, m_rlString.size(), (double)m_numSymbols / m_rlString.size());
    printf("Memory -- Markers: %zu Str: %zu Misc: %zu Total: %zu (%lf MB)\n", m_size, bwStr_size, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n", m_numSymbols, (double)total_size / m_numSymbols);
}

// Print the run length distribution of the BWT
void RLBWT::printRunLengths() const
{
    typedef std::map<size_t, size_t> DistMap;
    DistMap rlDist;

    size_t numRuns = getNumRuns();
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = m_rlString[i];
        size_t length = unit.getCount();
        rlDist[length]++;
    }
    
    printf("Run length distrubtion\n");
    printf("rl\tcount\tfrac\n");
    for(DistMap::iterator iter = rlDist.begin(); iter != rlDist.end(); ++iter)
    {
        printf("%zu\t%zu\t%lf\n", iter->first, iter->second, double(iter->second) / numRuns);
    }
    printf("Total runs: %zu\n", numRuns);
}

