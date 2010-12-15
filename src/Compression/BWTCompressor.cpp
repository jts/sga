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
                                 m_unitsWrote(0),
                                 m_pTempLargeMarkerWriter(NULL),
                                 m_pTempSmallMarkerWriter(NULL)
{

}

//
BWTCompressor::~BWTCompressor()
{
    std::cout << "BWTCompressor -- Wrote: " << m_symbolsWrote << " in " << m_unitsWrote << " bytes\n";
    std::cout << "BWTCompressor -- LargeMarkers placed: " << m_numLargeMarkersWrote << "\n";
    std::cout << "BWTCompressor -- SmallMarkers placed: " << m_numSmallMarkersWrote << "\n";
}

void BWTCompressor::initialize(const std::string& filename, size_t largeSampleRate, size_t smallSampleRate, size_t numSymbols)
{
    m_largeSampleRate = largeSampleRate; 
    m_smallSampleRate = smallSampleRate;
    m_expectedSymbols = numSymbols;

    m_numSmallMarkersWrote = 0;
    m_numLargeMarkersWrote = 0;

    // Open temporary files for the large and small markers
    assert(m_pTempLargeMarkerWriter == NULL && m_pTempSmallMarkerWriter == NULL);
    m_tempLargeFilename = filename + ".largemarkers.tmp";
    m_tempSmallFilename = filename + ".smallmarkers.tmp";
    m_pTempLargeMarkerWriter = createWriter(m_tempLargeFilename, std::ios::out | std::ios::binary);
    m_pTempSmallMarkerWriter = createWriter(m_tempSmallFilename, std::ios::out | std::ios::binary);

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

    // early exit if we do not have any data to write out
    if(m_symbolBuffer.empty())
        return;

    // Constrct new markers that are necessary
    buildMarkers();

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

    // Set the encoder index in the small marker and write it to the temp file
    m_prevSmallMarker.encoderIdx = encoderIdx;
    m_pTempSmallMarkerWriter->write(reinterpret_cast<const char*>(&m_prevSmallMarker), sizeof(SmallMarker));
    m_numSmallMarkersWrote += 1;

    // Perform the actual encoding
    ByteVector encodedBytes;
    size_t symbolsEncoded = StreamEncode::encodeStream(m_symbolBuffer, symbolEncoder, m_rlEncoder, encodedBytes);
    assert(symbolsEncoded == m_symbolBuffer.size());

    // Write the bytes out to the stream
    for(size_t i = 0; i < encodedBytes.size(); ++i)
        pWriter->write(reinterpret_cast<const char*>(&encodedBytes[i]), sizeof(uint8_t));

    m_unitsWrote += encodedBytes.size();
    m_symbolsWrote += symbolsEncoded;

    // Clear the buffer
    m_symbolBuffer.clear();
}

// Move the large markers from the temp files to pWriter
void BWTCompressor::writeLargeMarkers(std::ostream* pWriter)
{
    // Close the temporary file writer
    delete m_pTempLargeMarkerWriter;
    m_pTempLargeMarkerWriter = NULL;

    // Open a binary reader for the large markers
    std::istream* pReader = createReader(m_tempLargeFilename, std::ios::binary);
    
    // Write a header line to the outfile with the number of markers
    pWriter->write(reinterpret_cast<const char*>(&m_numLargeMarkersWrote), sizeof(m_numLargeMarkersWrote));

    // Copy the file byte by byte to pWriter
    binaryFileCopy(pReader, pWriter, m_numLargeMarkersWrote * sizeof(LargeMarker), sizeof(LargeMarker));
    delete pReader; 
    unlink(m_tempLargeFilename.c_str());
}

// Move the small markers from the temp files to pWriter
void BWTCompressor::writeSmallMarkers(std::ostream* pWriter)
{
    delete m_pTempSmallMarkerWriter;
    m_pTempSmallMarkerWriter = NULL;

    // Open a binary reader for the large markers
    std::istream* pReader = createReader(m_tempSmallFilename, std::ios::binary);
    
    // Write a header line to the outfile with the number of markers
    pWriter->write(reinterpret_cast<const char*>(&m_numSmallMarkersWrote), sizeof(m_numSmallMarkersWrote));

    // Copy the file byte by byte to pWriter
    binaryFileCopy(pReader, pWriter, m_numSmallMarkersWrote * sizeof(SmallMarker), sizeof(SmallMarker));
    delete pReader; 
    unlink(m_tempSmallFilename.c_str());
}


// Copy numBytes from pReader to pWriter in units of size unitSize
void BWTCompressor::binaryFileCopy(std::istream* pReader, std::ostream* pWriter, 
                                   size_t numBytes, size_t unitSize)
{
    char buffer[unitSize];
    assert(numBytes % unitSize == 0);
    size_t bytesCopied = 0;
    while(bytesCopied < numBytes)
    {
        assert(pReader->good());
        pReader->read(buffer, unitSize);
        pWriter->write(buffer, unitSize);
        bytesCopied += unitSize;
    }
}


// Construct markers for the sequence up to this point
void BWTCompressor::buildMarkers()
{
    // Large marker is placed if its been sampleRateLarge symbols since the last one
    size_t startingUnitIndex = m_unitsWrote;
    while((m_symbolsWrote / m_largeSampleRate) + 1 > m_numLargeMarkersWrote)
    {
        // Build a new large marker with the accumulated counts up to this point
        LargeMarker marker;
        marker.unitIndex = startingUnitIndex;
        marker.counts = m_runningAC;
        m_prevLargeMarker = marker;

        // Write the marker to the temp file
        m_pTempLargeMarkerWriter->write(reinterpret_cast<const char*>(&m_prevLargeMarker), sizeof(LargeMarker));
        m_numLargeMarkersWrote += 1;
    }

    // SmallMarkers are placed every segment. 
    AlphaCount16 smallAC;
    for(size_t j = 0; j < ALPHABET_SIZE; ++j)
    {
        size_t v = m_runningAC.getByIdx(j) - m_prevLargeMarker.counts.getByIdx(j);
        if(v > smallAC.getMaxValue())
        {
            std::cerr << "Error: Number of symbols exceeds the maximum value (" << v << " > " << smallAC.getMaxValue() << ")\n";
            std::cerr << "RunningAC: " << m_runningAC << "\n";
            std::cerr << "PrevAC: " << m_prevLargeMarker.counts << "\n";
            std::cerr << "SmallAC:" << smallAC << "\n";
            exit(EXIT_FAILURE);
        }
        smallAC.setByIdx(j, v);
    }
    
    // Construct the small marker
    SmallMarker smallMarker;
    smallMarker.unitCount = startingUnitIndex - m_prevLargeMarker.unitIndex;
    smallMarker.counts = smallAC;        
    m_prevSmallMarker = smallMarker;
}

AlphaCount64 BWTCompressor::getRunningCount() const
{
    return m_runningAC;
}

// 
size_t BWTCompressor::getNumBytesWrote() const
{
    return m_unitsWrote;
}

//
size_t BWTCompressor::getNumSymbolsWrote() const
{
    return m_symbolsWrote;
}
