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
    size_t maxSyms = 5;
    size_t minRun = 3;
    
    size_t numSymbolsWrote = 0;
    size_t numRLWrote = 0;
    size_t numHuffWrote = 0;

    typedef std::deque<char> CharDeque;
    CharDeque symbolBuffer;

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

        // Encode the current characters in the buffer
        bool encodeDone = false;
        while(!encodeDone)
        {
            if(symbolBuffer.size() < maxSyms)
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

            if(currRunLength >= minRun)
            {
                // Encode using RLE
                RLUnit currUnit(first);
                symbolBuffer.pop_front();
                size_t lenToEncode = std::min((size_t)RL_FULL_COUNT, currRunLength) - 1;
                for(size_t i = 0; i < lenToEncode; ++i)
                {
                    currUnit.incrementCount();
                    symbolBuffer.pop_front();
                }

                // Write out the unit
                out.push_back(currUnit.data);
                numSymbolsWrote += currUnit.getCount();
                numRLWrote += 1;
            }
            else
            {
                // Encode using Huffman
                size_t n = std::min(symbolBuffer.size(), maxSyms);
                if(n < maxSyms)
                {
                    // stop encoding if there aren't enough symbols to fill up
                    encodeDone = true;
                    break;
                }

                HuffUnit currUnit;
                for(size_t i = 0; i < n; ++i)
                {
                    char r = symbolBuffer.front();
                    currUnit.setChar(r, i);
                    symbolBuffer.pop_front();
                    numSymbolsWrote += 1;
                }

                // Push to the stream
                char* bytes = (char*)&currUnit.data;
                out.push_back(bytes[0]);
                out.push_back(bytes[1]);
                numHuffWrote += 1;
            }
        }
    }

    std::cout << "Wrote " << numSymbolsWrote << " symbols\n";
    std::cout << "Wrote " << numRLWrote << " runs\n";
    std::cout << "Wrote " << numHuffWrote << " huff\n";
    size_t bytes = numRLWrote + 2*numHuffWrote;
    std::cout << "Total: " << bytes << " bytes\n";
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
