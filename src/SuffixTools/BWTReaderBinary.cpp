//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLBWTReaderBinary - Read a run length encoded BWT file from disk
//
#include <algorithm>
#include "BWTReaderBinary.h"
#include "SBWT.h"
#include "RLBWT.h"
#include "Timer.h"
#include "HuffmanUtil.h"
#include "StreamEncoding.h"
#include "HuffmanForest.h"

static std::string hackname;
//
BWTReaderBinary::BWTReaderBinary(const std::string& filename) : m_stage(IOS_NONE), m_numUnitsOnDisk(0), m_numUnitsRead(0)
{
    m_pReader = createReader(filename, std::ios::binary);
    m_stage = IOS_HEADER;
    HuffmanTreeCodec<int> huffTree = HuffmanUtil::buildRLHuffmanTree();
    m_rlDecodeTable.initialize(huffTree);
    hackname = filename;

    // streaming mode state variables
    m_numUnitsRead = 0;
    m_numSymbolsTotal = 0;
    m_totalSmallMarkers = 0;
    m_readSmallMarkers = 0;
    m_smallMarkerOffset = 0;
    m_decompressedTotal = 0;
}

//
BWTReaderBinary::~BWTReaderBinary()
{
    delete m_pReader;
}

void BWTReaderBinary::read(RLBWT* pRLBWT)
{
    Timer readTimer("reader", false);
    BWFlag flag;

    readHeader(pRLBWT->m_numStrings, pRLBWT->m_numSymbols, flag);

    pRLBWT->setSampleRates(m_largeSampleRate, m_smallSampleRate);

    // Read the compressed data and markers from disk
    assert(m_numUnitsOnDisk > 0);
    readCompressedData(pRLBWT, m_numUnitsOnDisk);
    readMarkers(pRLBWT);
    readPredCounts(pRLBWT);

    //pRLBWT->printInfo();
    //pRLBWT->print();

    printf("Done reading BWT (%3.2lfs)\n", readTimer.getElapsedWallTime());
    pRLBWT->decodeToFile(hackname + ".decode");
}

void BWTReaderBinary::read(SBWT* /*pSBWT*/)
{
    // deprecated
    assert(false);
#if 0
    BWFlag flag;
    size_t numSymbols;
    size_t largeSampleRate = 0;
    size_t smallSampleRate = 0;

    readHeader(pSBWT->m_numStrings, numSymbols, largeSampleRate, smallSampleRate, flag);
    size_t numRead = 0;
    pSBWT->m_bwStr.resize(numSymbols);
    while(numRead < numSymbols)
    {
        char b = readBWChar();
        pSBWT->m_bwStr.set(numRead++, b);
    }

    assert(m_numRunsOnDisk > 0);
#endif
}

//
void BWTReaderBinary::readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    uint16_t magic_number;
    m_pReader->read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
    
    if(magic_number != RLBWT_FILE_MAGIC)
    {
        std::cerr << "BWT file is not properly formatted, aborting\n";
        exit(EXIT_FAILURE);
    }

    size_t largeMarkerOffset = 0;

    m_pReader->read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    m_pReader->read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    m_pReader->read(reinterpret_cast<char*>(&m_largeSampleRate), sizeof(m_largeSampleRate));
    m_pReader->read(reinterpret_cast<char*>(&m_smallSampleRate), sizeof(m_smallSampleRate));
    m_pReader->read(reinterpret_cast<char*>(&m_numUnitsOnDisk), sizeof(m_numUnitsOnDisk));
    m_pReader->read(reinterpret_cast<char*>(&largeMarkerOffset), sizeof(largeMarkerOffset));
    m_pReader->read(reinterpret_cast<char*>(&m_smallMarkerOffset), sizeof(m_smallMarkerOffset));
    m_pReader->read(reinterpret_cast<char*>(&flag), sizeof(flag));

    // Save the total number of symbols to read
    m_numSymbolsTotal = num_symbols;
    
    // When reading the BWT file one symbol at a time we need the small markers to decode each data segment.
    // We record the position in the file of the last small marker read and the number of small markers parsed so far
    std::streampos currPos = m_pReader->tellg();
    m_pReader->seekg(m_smallMarkerOffset);
    m_pReader->read(reinterpret_cast<char*>(&m_totalSmallMarkers), sizeof(m_totalSmallMarkers));
    m_smallMarkerOffset = m_pReader->tellg(); // m_smallMarkerOffset now points to the first small marker to read
    m_pReader->seekg(currPos);

    /*
    std::cout << "Read magic: " << magic_number << "\n";
    std::cout << "strings: " << num_strings << "\n";
    std::cout << "symbols: " << num_symbols << "\n";
    std::cout << "units: " << m_numUnitsOnDisk << "\n";
    std::cout << "large sample rate: " << m_largeSampleRate << "\n";
    std::cout << "small sample rate: " << m_smallSampleRate << "\n";
    std::cout << "large marker offset: " << largeMarkerOffset << "\n";
    std::cout << "small marker offset: " << m_smallMarkerOffset << "\n";
    std::cout << "num small markers: " << m_totalSmallMarkers << "\n";
    */
    m_stage = IOS_BWSTR;    
}

void BWTReaderBinary::readCompressedData(RLBWT* pBWT, size_t numUnits)
{
    pBWT->m_rlString.resize(numUnits);
    m_pReader->read(reinterpret_cast<char*>(&pBWT->m_rlString[0]), numUnits*sizeof(uint8_t));
}

void BWTReaderBinary::readMarkers(RLBWT* pBWT)
{
    // Read large markers
    size_t numLargeMarkers;
    m_pReader->read(reinterpret_cast<char*>(&numLargeMarkers), sizeof(numLargeMarkers));
    pBWT->m_largeMarkers.resize(numLargeMarkers);
    m_pReader->read(reinterpret_cast<char*>(&pBWT->m_largeMarkers[0]), sizeof(LargeMarker) * numLargeMarkers);

    size_t numSmallMarkers;
    m_pReader->read(reinterpret_cast<char*>(&numSmallMarkers), sizeof(numSmallMarkers));
    pBWT->m_smallMarkers.resize(numSmallMarkers);
    m_pReader->read(reinterpret_cast<char*>(&pBWT->m_smallMarkers[0]), sizeof(SmallMarker) * numSmallMarkers);
}

void BWTReaderBinary::readPredCounts(RLBWT* pBWT)
{
    // Read an alphacount from the stream
    m_pReader->read(reinterpret_cast<char*>(&pBWT->m_predCount), sizeof(pBWT->m_predCount));
}

// Read a single base from the BWStr
// The BWT is stored as a compressed sequence of bytes on the disk
// We decode entire units at once, by reading in a small marker and then at least
// m_smallSampleRate units from the disk file.
char BWTReaderBinary::readBWChar()
{
    assert(m_stage == IOS_BWSTR);

    if(m_decompressedBuffer.empty())
    {
        if(m_decompressedTotal == m_numSymbolsTotal)
        {
            assert(m_numUnitsRead == m_numUnitsOnDisk);
            assert(m_readSmallMarkers == m_totalSmallMarkers);
            return '\n';
        }
        else
        {
            // Read a block in and decompress it
            decompressBlock();
        }
    }

    assert(!m_decompressedBuffer.empty());
    char b = m_decompressedBuffer.front();
    m_decompressedBuffer.pop_front();
    return b;
}

// Read some data from the file and decompress it to the buffer
void BWTReaderBinary::decompressBlock()
{
    // Read in the next small marker
    SmallMarker smallMarker = readSmallMarker();

    // Read data in the buffer so that it contains at least m_smallSampleRate chars
    fillCompressedBuffer();

    // Perform the decompression
    _decompressBuffer(smallMarker);
}

// Read a single small marker from the stream
// We seek to the position in the stream of the next small marker, read it, then seek back
SmallMarker BWTReaderBinary::readSmallMarker()
{
    assert(m_readSmallMarkers < m_totalSmallMarkers);
    size_t currFilePosition = m_pReader->tellg();
    m_pReader->seekg(m_smallMarkerOffset);
    SmallMarker smallMarker;
    m_pReader->read(reinterpret_cast<char*>(&smallMarker), sizeof(smallMarker));
    m_smallMarkerOffset = m_pReader->tellg();
    m_pReader->seekg(currFilePosition);
    m_readSmallMarkers += 1;
    return smallMarker;
}

// Read from the compressed bwt until the compressed buffer has m_smallSampleRate character
void BWTReaderBinary::fillCompressedBuffer()
{
    size_t numUnitsToRead = std::min(m_smallSampleRate - m_compressedBuffer.size(), m_numUnitsOnDisk - m_numUnitsRead);
    //std::cout << "Reading " << numUnitsToRead << " from file\n";
    std::vector<uint8_t> readBuffer;
    readBuffer.resize(numUnitsToRead);
    m_pReader->read(reinterpret_cast<char*>(&readBuffer[0]), numUnitsToRead);
    
    m_numUnitsRead += numUnitsToRead;
    m_compressedBuffer.insert(m_compressedBuffer.end(), readBuffer.begin(), readBuffer.end());
}

// Decompress the current buffer using the smallMarker
void BWTReaderBinary::_decompressBuffer(const SmallMarker& smallMarker)
{
    size_t numSymbolsToDecompress = std::min(m_numSymbolsTotal - m_decompressedTotal, m_smallSampleRate);
    StreamEncode::CharDequeDecode cdd(m_decompressedBuffer);
    size_t numBitsRead = 0;
    size_t numSymbols = StreamEncode::decodeStream(HuffmanForest::Instance().getDecoder(smallMarker.encoderIdx), 
                                                   m_rlDecodeTable, &m_compressedBuffer[0], (&m_compressedBuffer.back()) + 1, 
                                                   numSymbolsToDecompress, numBitsRead, cdd);

    assert(numSymbols == numSymbolsToDecompress);

    // Compute the number of bytes that were processed from the stream
    size_t numBytesRead = (numBitsRead % 8 == 0) ? numBitsRead / 8 : (numBitsRead / 8) + 1;
    assert(numBytesRead <= m_compressedBuffer.size());

    // Discard numBytesRead from the compressed buffer
    m_compressedBuffer.erase(m_compressedBuffer.begin(), m_compressedBuffer.begin() + numBytesRead);
    
    m_decompressedTotal += numSymbols;
}

