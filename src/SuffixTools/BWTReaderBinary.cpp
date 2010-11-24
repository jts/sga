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
    readRuns(pRLBWT->m_rlString, m_numRunsOnDisk);

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

void BWTReaderBinary::readRuns(RLRawData& out, size_t numRuns)
{
    //out.resize(numRuns);
    //m_pReader->read(reinterpret_cast<char*>(&out[0]), numRuns*sizeof(RLUnit));
    (void)numRuns;   
    bool readDone = false;
    //size_t minRun = 3;
    
    size_t numSymbolsWrote = 0;
    size_t numRLWrote = 0;
    size_t numHuffWrote = 0;

    size_t numSymbolsEncoded = 0;
    size_t totalBitsUsed = 0;
    size_t symbolBitsUsed = 0;
    size_t rlBitsUsed = 0;
    
    std::map<int, int> runBits;
    std::map<char, int> symbolBits;

    typedef std::deque<char> CharDeque;
    CharDeque symbolBuffer;
    
    HuffmanRunLengthEncoder runEncoder;
    Huffman::buildRLHuffmanTree(runEncoder);

    // Read at least 256 symbols from the stream as long as symbols remain to be read
    // We stop the read once 256 symbols are parsed and we are not in the middle of a run
    // Then we encode the string in a set of units and write them out
    size_t minBufferSize = 256;
    while(!readDone)
    {
        // Read symbols into the buffer
        bool inRun = true;
        while(symbolBuffer.size() < minBufferSize || inRun)
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

        // Count symbols in the buffer
        AlphaCount64 ac;
        for(size_t i = 0; i < symbolBuffer.size(); ++i)
        {
            ac.increment(symbolBuffer[i]);
        }

        // Get a sorting of the characters
        HuffmanSymbolEncoder symbolEncoder;
        HuffmanSymbolDecoder symbolDecoder;

        Huffman::buildSymbolHuffman(ac, symbolEncoder, symbolDecoder);

        // Encode the current characters in the buffer
        bool encodeDone = false;
        while(!encodeDone)
        {
            if(symbolBuffer.empty())
            {
                encodeDone = true;
                break;
            }

            // Count the length of the run at the front
            assert(!symbolBuffer.empty());
            char first = symbolBuffer.front();
            size_t currRunLength = 1;
            for(size_t i = 1; i < symbolBuffer.size(); ++i)
            {
                if(symbolBuffer[i] == first)
                    currRunLength += 1;
                else
                    break;
            }

            size_t encodingRunLength = Huffman::getEncodingLength(runEncoder, currRunLength);
            assert(encodingRunLength <= currRunLength);

            HuffmanEncodePair symPair = symbolEncoder[first];
            assert(runEncoder.find(encodingRunLength) != runEncoder.end());
            HuffmanEncodePair runPair = runEncoder[encodingRunLength];
            totalBitsUsed += symPair.bits;
            totalBitsUsed += runPair.bits;
            
            symbolBitsUsed += symPair.bits;
            rlBitsUsed += runPair.bits;
            numSymbolsWrote += encodingRunLength;
            numSymbolsEncoded += 1;
            runBits[encodingRunLength] += runPair.bits;
            symbolBits[first] += symPair.bits;

            // Remove the symbols encoded from the buffer
            for(size_t i = 0; i < encodingRunLength; ++i)
                symbolBuffer.pop_front();

            /*
            // Encode using Huffman
            uint8_t data = 0;
            size_t totalBits = 8;
            size_t totalBytes = totalBits / 8;
            size_t targetBit = 0;
            size_t numBitsRemaining = totalBits;
            size_t numSymbolsTaken = 0;
            while(numSymbolsTaken < symbolBuffer.size() && numBitsRemaining >= 4)
            {
                char r = symbolBuffer[numSymbolsTaken];
                assert(encoder.find(r) != encoder.end());
                HuffmanEncodePair pair = encoder[r];
                size_t code = pair.code;
                size_t codeBits = pair.bits;

                size_t currBit = totalBits - codeBits;
                size_t shift = currBit - targetBit;
                code <<= shift;
                data |= code;
                targetBit += codeBits;
                numBitsRemaining -= codeBits;
                numSymbolsTaken += 1;
            }

            // Ensure an entire unit was encoded
            if(numBitsRemaining < 4)
            {
                // Push to the stream
                char* bytes = (char*)&data;
                for(size_t i = 0; i < totalBytes; ++i)
                    out.push_back(bytes[i]);

                for(size_t i = 0; i < numSymbolsTaken; ++i)
                    symbolBuffer.pop_front();
                numSymbolsWrote += numSymbolsTaken;
                numHuffWrote += totalBytes;
            }
            else
            {
                encodeDone = true;
            }
            */
        }
    }

    for(std::map<int,int>::iterator iter = runBits.begin(); iter != runBits.end(); ++iter)
    {
        std::cout << iter->first << "\t" << iter->second / 8 << "\n";
    }
    std::cout << "Wrote " << numSymbolsWrote << " symbols\n";
    std::cout << "Wrote " << numRLWrote << " runs\n";
    std::cout << "Wrote " << numHuffWrote << " huff\n";
    //size_t bytes = out.size();
    (void)out;
    size_t bytes = totalBitsUsed / 8;
    std::cout << "Total: " << bytes << " bytes\n";
    std::cout << "BITS: " << totalBitsUsed << "\n";
    std::cout << "SYMBITS: " << symbolBitsUsed << " for " << numSymbolsEncoded << " (" << (double)symbolBitsUsed / numSymbolsEncoded << ")\n";
    std::cout << "RLEBITS: " << rlBitsUsed << " for " << numSymbolsEncoded << " (" << (double)rlBitsUsed / numSymbolsEncoded << ")\n";
    std::cout << (double)numSymbolsWrote / bytes << " symbols/byte\n";


    exit(1);
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
