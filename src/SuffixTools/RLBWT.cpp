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
#include "HuffmanForest.h"
#include "StreamEncoding.h"
#include <istream>
#include <queue>
#include <inttypes.h>

// macros
#define OCC(c,i) m_occurrence.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

// Parse a BWT from a file
RLBWT::RLBWT(const std::string& filename, int sampleRate) : m_numStrings(0), 
                                                            m_numSymbols(0), 
                                                            m_largeSampleRate(DEFAULT_SAMPLE_RATE_LARGE),
                                                            m_smallSampleRate(sampleRate)
{
    m_rlHuffman = HuffmanUtil::buildRLHuffmanTree();
    m_rlDecodeTable.initialize(m_rlHuffman);

    IBWTReader* pReader = BWTReader::createReader(filename);
    pReader->read(this);

    //initializeFMIndex();
    delete pReader;
}

//
void RLBWT::append(char b)
{
    assert(false && "deprecated");
    bool increment = false;
    if(!m_rlString.empty())
    {
        RLUnit* pLastUnit = (RLUnit*)&m_rlString.back();
        if(pLastUnit->getChar() == b && !pLastUnit->isFull())
        {
            pLastUnit->incrementCount();
            increment = true;
        }
    }

    if(!increment)
    {
        // Must add a new unit to the string
        RLUnit unit(b);
        m_rlString.push_back(unit.data);
    }
    ++m_numSymbols;
}

void RLBWT::setSampleRates(size_t largeSampleRate, size_t smallSampleRate)
{
    m_smallSampleRate = smallSampleRate;
    m_largeSampleRate = largeSampleRate;

    m_smallShiftValue = Occurrence::calculateShiftValue(m_smallSampleRate);
    m_largeShiftValue = Occurrence::calculateShiftValue(m_largeSampleRate);
}

// get the number of markers required to cover the n symbols at sample rate of d
size_t RLBWT::getNumRequiredMarkers(size_t n, size_t d) const
{
    // we place a marker at the beginning (with no accumulated counts), every m_sampleRate
    // bases and one at the very end (with the total counts)
    size_t num_markers = (n % d == 0) ? (n / d) + 1 : (n / d) + 2;
    return num_markers;
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
    size_t small_m_size = m_smallMarkers.size() * sizeof(SmallMarker);
    size_t large_m_size = m_largeMarkers.size() * sizeof(LargeMarker);
    size_t total_marker_size = small_m_size + large_m_size;

    size_t bwStr_size = m_rlString.size();
    size_t other_size = sizeof(*this);
    size_t total_size = total_marker_size + bwStr_size + other_size;

    double mb = (double)(1024 * 1024);
    double total_mb = total_size / mb;
    
    printf("\nRLBWT info:\n");
    printf("Large Sample rate: %zu\n", m_largeSampleRate);
    printf("Small Sample rate: %zu\n", m_smallSampleRate);
    printf("Contains %zu symbols in %zu bytes (%1.4lf symbols per bytes)\n", m_numSymbols, m_rlString.size(), (double)m_numSymbols / m_rlString.size());
    printf("Marker Memory -- Small Markers: %zu (%.1lf MB) Large Markers: %zu (%.1lf MB)\n", small_m_size, small_m_size / mb, large_m_size, large_m_size / mb);
    printf("Total Memory -- Markers: %zu (%.1lf MB) Str: %zu (%.1lf MB) Misc: %zu Total: %zu (%lf MB)\n", total_marker_size, total_marker_size / mb, bwStr_size, bwStr_size / mb, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n\n", m_numSymbols, (double)total_size / m_numSymbols);
}

// Print the run length distribution of the BWT
void RLBWT::printRunLengths() const
{
    typedef std::map<size_t, size_t> DistMap;
    DistMap rlDist;

    char prevSym = '\0';
    size_t prevRunLen = 0;
    size_t currLen = 0;
    size_t adjacentSingletons = 0;
    size_t numRuns = getNumRuns();
    size_t totalRuns = 0;
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = m_rlString[i];
        size_t length = unit.getCount();
        if(unit.getChar() == prevSym)
        {
            currLen += unit.getCount();
        }
        else
        {
            if(prevSym != '\0')
            {
                if(currLen >= 200)
                    rlDist[200]++;
                else
                    rlDist[currLen]++;
                totalRuns++;
            }
            currLen = length;
            prevSym = unit.getChar();
        }

        if(length == 1 && prevRunLen == 1)
        {
            adjacentSingletons += 1;
            prevRunLen = 0;
        }
        prevRunLen = length;
    }
    
    printf("Run length distrubtion\n");
    printf("rl\tcount\tfrac\n");
    double cumulative_mb = 0.0f;
    for(DistMap::iterator iter = rlDist.begin(); iter != rlDist.end(); ++iter)
    {
        cumulative_mb += ((double)iter->second / (1024 * 1024));
        printf("%zu\t%zu\t%lf\t%lf\n", iter->first, iter->second, double(iter->second) / totalRuns, cumulative_mb);
    }
    printf("Total runs: %zu\n", totalRuns);
    printf("Number of adjacent singleton runs: %zu\n", adjacentSingletons);
    printf("Minimal runs: %zu\n", totalRuns - adjacentSingletons / 2);
}

void RLBWT::decodeToFile(const std::string& filename)
{
    std::ofstream outFile(filename.c_str());
    size_t numSymbolsDecoded = 0;
    AlphaCount64 runningAC;

    // Iterate over the small markers, decoding each segment that the marker represents
    for(size_t i = 0; i < m_smallMarkers.size(); ++i)
    {
        size_t decodeIdx = 0;
        LargeMarker currLargeMarker = getInterpolatedMarker(i, decodeIdx);

        // Calculate the number of symbols to decode in the block
        size_t numToDecode = 0;
        if(i == (m_smallMarkers.size() - 1))
        {
             // last marker, decode the remainder of the string
             numToDecode = m_numSymbols - numSymbolsDecoded;
        }
        else
        {
            size_t nlmDecode = 0;
            LargeMarker nextLargeMarker = getInterpolatedMarker(i+1, nlmDecode);
            //std::cout << "NP: " << nextLargeMarker.getActualPosition() << " CP: " << currLargeMarker.getActualPosition() << "\n";
            numToDecode = nextLargeMarker.getActualPosition() - currLargeMarker.getActualPosition();
        }

        std::string outStr;
        const CharPackedTableDecoder& symDecoder = HuffmanForest::Instance().getDecoder(decodeIdx);
        size_t startingUnitIndex = currLargeMarker.unitIndex;
        StreamEncode::StringDecode sd(outStr);
        size_t numBitsRead = 0;
        size_t numDecoded = StreamEncode::decodeStream(symDecoder, m_rlDecodeTable, &m_rlString[startingUnitIndex], &m_rlString.back(), numToDecode, numBitsRead, sd);
        assert(numDecoded >= numToDecode);
        
        for(size_t j = 0; j < outStr.size(); ++j)
        {
            runningAC.increment(outStr[j]);
            size_t true_pos = numSymbolsDecoded + j;
            if(outStr[j] != getChar(true_pos))
            {
                std::cout << "String mismatch at j: " << j << ", truepos:" << true_pos << " " << outStr[j] << " != " << getChar(true_pos) << "\n";
            }
            AlphaCount64 calcAC = getFullOcc(true_pos);
            if(calcAC != runningAC)
            {
                std::cout << "Fail for idx: " << true_pos << "\n";
                std::cout << "CalcAC: " << calcAC << "\n";
                std::cout << "RunnAC: " << runningAC << "\n";
                exit(1);
            }
        }
        numSymbolsDecoded += numDecoded;
        outFile << i << " " << numBitsRead << " " << numDecoded << " " << outStr << "\n";
        //outFile << outStr << "\n";
    }
    assert(numSymbolsDecoded == m_numSymbols);
}

