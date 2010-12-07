//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTCompressor - Compress a BWT to disk
//
#include "BWTCompressor.h"
#include "HuffmanUtil.h"
#include "HuffmanForest.h"
#include "StreamEncoding.h"

BWTCompressor::BWTCompressor() : m_largeSampleRate(0),
                                 m_smallSampleRate(0), 
                                 m_expectedSymbols(0), 
                                 m_symbolsWrote(0), 
                                 m_unitsWrote(0)

{

}

//
BWTCompressor::~BWTCompressor()
{
    std::cout << "BWTCompressor -- Wrote: " << m_symbolsWrote << " in " << m_unitsWrote << " bytes\n";
    std::cout << "BWTCompressor -- LargeMarkers placed: " << m_largeMarkers.size() << "\n";
    std::cout << "BWTCompressor -- SmallMarkers placed: " << m_smallMarkers.size() << "\n";
}

void BWTCompressor::initialize(size_t largeSampleRate, size_t smallSampleRate, size_t numSymbols)
{
    m_largeSampleRate = largeSampleRate; 
    m_smallSampleRate = smallSampleRate;
    m_expectedSymbols = numSymbols;

    WARN_ONCE("BWTCompressor:: Resize marker arrays");

    // Initialize the run length encoder
    m_rlEncoder = HuffmanUtil::buildRLHuffmanTree();

}

// Add a symbol to the buffer and flush it to disk if necessary
void BWTCompressor::writeSymbol(char symbol, std::ostream* pWriter)
{
    assert(m_smallSampleRate > 0);
    assert(m_largeSampleRate > 0);
    assert(m_expectedSymbols > 0);

    m_symbolBuffer.push_back(symbol);
    if(m_symbolBuffer.size() == m_smallSampleRate)
        flush(pWriter);
}

// Flush the current buffer to disk
void BWTCompressor::flush(std::ostream* pWriter)
{
    assert(m_symbolBuffer.size() <= m_smallSampleRate);

    // Add any new markers that are necessary
    SmallMarker& smallMarker = appendMarkers();

    // Count the number of runs of a given symbol in the buffer
    // Also, update the running count
    CharIntMap symbolCounts;
    char prevChar = '\0';
    for(size_t i = 0; i < m_symbolBuffer.size(); ++i)
    {
        char c = m_symbolBuffer[i];
        // skip runs to get an accurate count for the huffman
        if(c != prevChar)
            symbolCounts[c]++;
        prevChar = c;
        m_runningAC.increment(c);
    }

    // Based on the symbol counts, choose a huffman tree from the symbols
    size_t encoderIdx = 0;
    const HuffmanTreeCodec<char>& symbolEncoder = HuffmanForest::Instance().getBestEncoder(symbolCounts, encoderIdx);

    // Set the encoder index in the small marker
    smallMarker.encoderIdx = encoderIdx;

    // Perform the actual encoding
    ByteVector encodedBytes;
    size_t symbolsEncoded = StreamEncode::encodeStream(m_symbolBuffer, symbolEncoder, m_rlEncoder, encodedBytes);
    assert(symbolsEncoded == m_symbolBuffer.size());

    // Write the bytes out to the stream
    m_unitsWrote += encodedBytes.size();
    m_symbolsWrote += symbolsEncoded;
    (void)pWriter;

    // Clear the buffer
    m_symbolBuffer.clear();

}

// Add markers
SmallMarker& BWTCompressor::appendMarkers()
{
    // Large marker is placed if its been sampleRateLarge symbols since the last one
    size_t startingUnitIndex = m_unitsWrote;
    while((m_symbolsWrote / m_largeSampleRate) + 1 > m_largeMarkers.size())
    {
        // Place a new large marker with the accumulated counts up to this point
        LargeMarker marker;
        marker.unitIndex = startingUnitIndex;
        marker.counts = m_runningAC;
        m_largeMarkers.push_back(marker);
        //std::cout << "Placed large marker at " << numSymbolsWrote << " idx: " << startingUnitIndex << " MI: " << numSymbolsWrote / sampleRateLarge << "\n";
    }

    // SmallMarkers are placed every segment. 
    // Calculate the large marker that this small marker is relative to
    size_t smallMarkerPosition = m_smallMarkers.size() * m_smallSampleRate;
    size_t largeMarkerIdx = smallMarkerPosition / m_largeSampleRate;
    assert(largeMarkerIdx == m_largeMarkers.size() - 1);
    LargeMarker& prevLargeMarker = m_largeMarkers[largeMarkerIdx];
   
    //std::cout << "Placing small marker " << pBWT->m_smallMarkers.size() << " relative to LI: " << largeMarkerIdx << "\n";
    //std::cout << "Small marker pos: " << smallMarkerPosition << " largeMarkerPos: " << prevLargeMarker.getActualPosition() << "\n";
    AlphaCount16 smallAC;
    for(size_t j = 0; j < ALPHABET_SIZE; ++j)
    {
        size_t v = m_runningAC.getByIdx(j) - prevLargeMarker.counts.getByIdx(j);
        if(v > smallAC.getMaxValue())
        {
            std::cerr << "Error: Number of symbols exceeds the maximum value (" << v << " > " << smallAC.getMaxValue() << ")\n";
            std::cerr << "RunningAC: " << m_runningAC << "\n";
            std::cerr << "PrevAC: " << prevLargeMarker.counts << "\n";
            std::cerr << "SmallAC:" << smallAC << "\n";
            exit(EXIT_FAILURE);
        }
        smallAC.setByIdx(j, v);
    }
    
    // Set the small marker
    SmallMarker smallMarker;
    smallMarker.unitCount = startingUnitIndex - prevLargeMarker.unitIndex;
    smallMarker.counts = smallAC;        
    m_smallMarkers.push_back(smallMarker);
    return m_smallMarkers.back();
}

// 
size_t BWTCompressor::getNumBytesWrote() const
{
    return m_unitsWrote;
}

//
const LargeMarkerVector& BWTCompressor::getLargeMarkerVector() const
{
    return m_largeMarkers;
}

//
const SmallMarkerVector& BWTCompressor::getSmallMarkerVector() const
{
    return m_smallMarkers;
}
