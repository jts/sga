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
#include "Huffman.h"
#include "StreamEncoding.h"
#include "HuffmanForest.h"

//
BWTReaderBinary::BWTReaderBinary(const std::string& filename) : m_stage(IOS_NONE), m_numRunsOnDisk(0), m_numRunsRead(0)
{
    m_pReader = createReader(filename, std::ios::binary);
    m_stage = IOS_HEADER;
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

    assert(m_numRunsOnDisk > 0);
    readRuns(pRLBWT, pRLBWT->m_rlString, m_numRunsOnDisk);

    //pRLBWT->printInfo();
    //pRLBWT->print();
    printf("Done reading BWT (%3.2lfs)\n", readTimer.getElapsedWallTime());
}

void BWTReaderBinary::read(SBWT* pSBWT)
{
    BWFlag flag;
    size_t numSymbols;
    readHeader(pSBWT->m_numStrings, numSymbols, flag);
    size_t numRead = 0;
    pSBWT->m_bwStr.resize(numSymbols);
    while(numRead < numSymbols)
    {
        char b = readBWChar();
        pSBWT->m_bwStr.set(numRead++, b);
    }

    assert(m_numRunsOnDisk > 0);
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

    m_pReader->read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    m_pReader->read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    m_pReader->read(reinterpret_cast<char*>(&m_numRunsOnDisk), sizeof(m_numRunsOnDisk));
    m_pReader->read(reinterpret_cast<char*>(&flag), sizeof(flag));
    
    //std::cout << "Read magic: " << magic_number << "\n";
    //std::cout << "strings:" << num_strings << "\n";
    //std::cout << "symbols: " << num_symbols << "\n";
    //std::cout << "runs: " << m_numRunsOnDisk << "\n";
    
    m_stage = IOS_BWSTR;    
}

void BWTReaderBinary::readRuns(RLBWT* pBWT, RLRawData& out, size_t numRuns)
{
    //out.resize(numRuns);
    //m_pReader->read(reinterpret_cast<char*>(&out[0]), numRuns*sizeof(RLUnit));
    (void)numRuns;   
    bool readDone = false;
    
    size_t sampleRateLarge = RLBWT::DEFAULT_SAMPLE_RATE_LARGE;
    size_t sampleRateSmall = RLBWT::DEFAULT_SAMPLE_RATE_SMALL;

    size_t numSymbolsWrote = 0;
    size_t numBytesUsed = 0;
    size_t numSymbolsEncoded = 0;
    AlphaCount64 runningAC;

    typedef std::map<char, int> CharIntMap;
    CharIntMap symbolCounts;
    CharDeque symbolBuffer;
    
    HuffmanTreeCodec<int> rlEncoder = Huffman::buildRLHuffmanTree();
    std::cout << "Huffman RL encoder needs at most " << rlEncoder.getMaxBits() << " bits\n";

    // Read symbols from the stream into the buffer as long as symbols remain to be read
    // We do not want to break up runs so the buffer size might be slightly larger than the target
    size_t minBufferSize = sampleRateSmall;
    while(!readDone)
    {
        // Read symbols into the buffer
        bool inRun = true;
        while(symbolBuffer.size() < minBufferSize)
        {
            char b = readBWChar();
            if(b == '\n')
            {
                readDone = true;
                break;
            }

            // If this symbol is the same as the current end symbol, set the flag
            // that it is a run so we continue to read
            if(symbolBuffer.empty() || symbolBuffer.back() == b)
                inRun = true;
            else
                inRun = false;
            symbolBuffer.push_back(b);
        }
        
        //if(symbolBuffer.empty() && readDone)
        //    break;

        //
        // Write markers at the beginning of this data segment
        //

        // Large marker is placed if its been sampleRateLarge symbols since the last one
        size_t startingUnitIndex = pBWT->m_rlString.size();
        while((numSymbolsWrote / sampleRateLarge) + 1 > pBWT->m_largeMarkers.size())
        {
            // Place a new large marker with the accumulated counts up to this point
            LargeMarker marker;
            marker.unitIndex = startingUnitIndex;
            marker.counts = runningAC;
            pBWT->m_largeMarkers.push_back(marker);
            //std::cout << "Placed large marker at " << numSymbolsWrote << " idx: " << startingUnitIndex << " MI: " << numSymbolsWrote / sampleRateLarge << "\n";
        }

        // SmallMarkers are placed every segment. 
        // Calculate the large marker that this small marker is relative to
        size_t smallMarkerPosition = pBWT->m_smallMarkers.size() * sampleRateSmall;
        size_t largeMarkerIdx = smallMarkerPosition / sampleRateLarge;
        assert(largeMarkerIdx < pBWT->m_largeMarkers.size());
        LargeMarker& prevLargeMarker = pBWT->m_largeMarkers[largeMarkerIdx];

        //std::cout << "Placing small marker " << pBWT->m_smallMarkers.size() << " relative to LI: " << largeMarkerIdx << "\n";
        //std::cout << "Small marker pos: " << smallMarkerPosition << " largeMarkerPos: " << prevLargeMarker.getActualPosition() << "\n";
        AlphaCount16 smallAC;
        for(size_t j = 0; j < ALPHABET_SIZE; ++j)
        {
            size_t v = runningAC.getByIdx(j) - prevLargeMarker.counts.getByIdx(j);
            if(v > smallAC.getMaxValue())
            {
                std::cerr << "Error: Number of symbols exceeds the maximum value (" << v << " > " << smallAC.getMaxValue() << ")\n";
                std::cerr << "RunningAC: " << runningAC << "\n";
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
        pBWT->m_smallMarkers.push_back(smallMarker);

        // Count symbols in the buffer
        symbolCounts.clear();
        char prevChar = '\0';
        for(size_t i = 0; i < symbolBuffer.size(); ++i)
        {
            // skip runs to get an accurate count for the huffman
            if(symbolBuffer[i] != prevChar)
            {
                symbolCounts[symbolBuffer[i]]++;
            }
            prevChar = symbolBuffer[i];
            runningAC.increment(symbolBuffer[i]);
        }

        // Construct the huffman tree for the symbols based on the counts
        size_t encoderIdx = 0;
        const HuffmanTreeCodec<char>& symbolEncoder = HuffmanForest::Instance().getBestEncoder(symbolCounts, encoderIdx);

        // Set the encoder index in the last small marker
        pBWT->m_smallMarkers.back().encoderIdx = encoderIdx;

        // Perform the actual encoding
        EncodedArray encodedBytes;
        size_t symbolsEncoded = StreamEncode::encodeStream(symbolBuffer, symbolEncoder, rlEncoder, encodedBytes);
        assert(symbolsEncoded == symbolBuffer.size());

        // Copy the encoded bytes
        for(size_t i = 0; i < encodedBytes.size(); ++i)
        {
            pBWT->m_rlString.push_back(encodedBytes[i]);
        }
        numSymbolsWrote += symbolsEncoded;
        numBytesUsed += encodedBytes.size();
        numSymbolsEncoded += symbolsEncoded;
        
        /*
        std::string testOut;
        StreamEncode::StringDecode stringDecoder(testOut);
        StreamEncode::decodeStream(symbolEncoder, rlEncoder, &encodedBytes[0], symbolsEncoded, stringDecoder);
        
        // Validate the decoding
        bool failedValidate = false;
        if(testOut.size() != symbolBuffer.size())
        {
            failedValidate = true;
        }
        else
        {
            for(size_t i = 0; i < testOut.size(); ++i)
            {
                if(testOut[i] != symbolBuffer[i])
                    failedValidate = true;
            }
        }

        if(failedValidate)
        {
            std::cout << "Failed at " << numSymbolsEncoded << "\n";
            std::cout << "Failure to decode buffer of size " << symbolBuffer.size() << "\n";
            std::cout << "Decode size: " << testOut.size() << "\n";
            std::cout << "Input:   ";

            for(size_t i = 0; i < symbolBuffer.size(); ++i)
            {
                std::cout << symbolBuffer[i];
            }
            std::cout << "\nDecoded: " << testOut << "\n";
            assert(false);
        }
        */
        symbolBuffer.clear();
    }
    pBWT->initializeFMIndex(runningAC);
    //pBWT->decodeToFile("test.txt");
    (void)out;
    std::cout << "Wrote " << numSymbolsWrote << " symbols\n";
    //size_t bytes = out.size();
    (void)out;
    size_t bytes = numBytesUsed;
    size_t bits = bytes * 8;
    double bitsPerSym = ((double)bits / numSymbolsEncoded);
    double symPerByte = ((double)numSymbolsEncoded / bytes);

    std::cout << "Total: " << bytes << " bytes\n";
    std::cout << bitsPerSym << " bits/sym\n"; 
    std::cout << symPerByte << " sym/byte\n"; 
    std::cout << "RLStrLen: " << pBWT->m_rlString.size() << "\n";
//    std::cout << "SYMBITS: " << symbolBitsUsed << " for " << numSymbolsEncoded << " (" << (double)symbolBitsUsed / numSymbolsEncoded << ")\n";
//    std::cout << "RLEBITS: " << rlBitsUsed << " for " << numSymbolsEncoded << " (" << (double)rlBitsUsed / numSymbolsEncoded << ")\n";
    pBWT->initializeFMIndex(runningAC);
}

// Read a single base from the BWStr
// The BWT is stored as runs on disk, so this class keeps
// an internal buffer of a single run and emits characters from this buffer
// and performs reads as necessary. If all the runs have been read, emit
// a newline character to signal the end of the BWT
char BWTReaderBinary::readBWChar()
{
    assert(m_stage == IOS_BWSTR);

    if(m_currRun.isEmpty())
    {
        // All runs have been read and emitted, return the end marker
        if(m_numRunsRead == m_numRunsOnDisk)
            return '\n';
 
        // Read one run from disk
        m_pReader->read(reinterpret_cast<char*>(&m_currRun), sizeof(RLUnit));
        ++m_numRunsRead;
    }

    // Decrement the current run and emit its symbol
    m_currRun.decrementCount();
    return m_currRun.getChar();
}
