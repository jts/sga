//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// StreamEncoding -- Encode a stream of symbols
// using a huffman encoder. 
//
#ifndef STREAMENCODING_H
#define STREAMENCODING_H

#include "Occurrence.h"
#include "HuffmanUtil.h"
#include "PackedTableDecoder.h"
#include "RLE.h"

//#define DEBUG_ENCODING 1
#define BITS_PER_BYTE 8

typedef std::pair<char, int> RLEPair;
typedef std::vector<RLEPair> RLEPairVector;
typedef std::deque<char> CharDeque;
typedef std::vector<uint8_t> ByteVector;

namespace StreamEncode
{

    // Decode functors for the generic decoding function
    struct AlphaCountDecode
    {
        AlphaCountDecode(AlphaCount64& target) : m_target(target) {}
        inline void operator()(int rank, int rl)
        {
            m_target.addByIdx(rank, rl);
        }
        AlphaCount64& m_target;
    };

    struct StringDecode
    {
        StringDecode(std::string& target) : m_target(target) {}
        inline void operator()(int rank, int rl)
        {
            m_target.append(rl, BWT_ALPHABET::getChar(rank));
        }
        std::string& m_target;
    };

    struct CharDequeDecode
    {
        CharDequeDecode(CharDeque& target) : m_target(target) {}
        inline void operator()(int rank, int rl)
        {
            for(int i = 0; i < rl; ++i)
                m_target.push_back(BWT_ALPHABET::getChar(rank));
        }
        CharDeque& m_target;
    };
    
    struct BaseCountDecode
    {
        BaseCountDecode(char targetBase, size_t& targetCount) : m_targetBase(targetBase), 
                                                                    m_targetCount(targetCount) {}
        inline void operator()(int rank, int rl)
        {
            if(BWT_ALPHABET::getChar(rank) == m_targetBase)
                m_targetCount += rl;
        }
        char m_targetBase;
        size_t& m_targetCount;
    };   

    // Decoder which returns the last base added. This is used to extract a particular character from the stream
    struct SingleBaseDecode
    {
        SingleBaseDecode(char& base) : m_base(base) {}
        inline void operator()(int rank, int /*rl*/)
        {
            m_base = BWT_ALPHABET::getChar(rank);
        }
        char& m_base;
    };

    //
    inline void printEncoding(const ByteVector& output)
    {
        std::cout << "Encoding: ";
        for(size_t i = 0; i < output.size(); ++i)
        {
            std::cout << int2Binary(output[i], 8) << " ";
        }
        std::cout << "\n";
    }

    // Write the code into the output stream starting at currBit
    inline size_t _writeCode(int code, int bits, size_t currBit, ByteVector& output)
    {
#ifdef DEBUG_ENCODING
        printEncoding(output);
        std::cout << "Writing the code " << int2Binary(code, bits) << " at bit " << currBit << "\n";
#endif

        int bitsRemaining = bits;
        int bitOffset = 0;
        //std::cout << "WRITE CODE: " << int2Binary(code, bits) << "\n";
        while(bitsRemaining > 0)
        {
            // Calculate position to start the write
            int byte = currBit / BITS_PER_BYTE;
            int bitIdx = MOD_POWER_2(currBit, BITS_PER_BYTE);
            int bitsToWrite = std::min((BITS_PER_BYTE - bitIdx), bitsRemaining);

            // Calculate the shift values and masks to apply
            int currPos = (BITS_PER_BYTE - bits + bitOffset);
            
            // Mask off the bits we want
            int mask = ((1 << bitsToWrite) - 1) << (bits - (bitsToWrite + bitOffset));
            int inCode = code & mask;

#ifdef DEBUG_ENCODING
            std::cout << "Mask: " << int2Binary(mask, bits) << "\n";
            std::cout << "Masked: " << int2Binary(inCode, bits) << "\n";
#endif

            // Shift the code into position
            if(currPos < bitIdx)
                inCode >>= (bitIdx - currPos);
            else if(currPos > bitIdx)
                inCode <<= (currPos - bitIdx);

#ifdef DEBUG_ENCODING
            std::cout << "Shifted: " << int2Binary(inCode,8) << "\n";
#endif
            // set the value with an OR
            output[byte] |= inCode;

            bitsRemaining -= bitsToWrite;
            bitOffset += bitsToWrite;
            currBit += bitsToWrite;
        }
        return bits;
    }

    // Read maxBits from the array starting at currBits and write the value to outCode
    // Returns the number of bits read
    inline void _readCode(const uint16_t currBit, const uint16_t baseShift, const uint16_t mask, const uint16_t input, uint16_t& outCode)
    {
#ifdef DEBUG_ENCODING
        printEncoding(input);
        std::cout << "Reading " << maxBits << " from array starting at " << currBit << "\n";
        std::cout << "Mask " << int2Binary(mask, maxBits) << "\n";
#endif      
        outCode = (input >> (baseShift - currBit)) & mask;
    }
    
    template<typename Functor>
    size_t decodeStream(const unsigned char* pInput, const unsigned char* /*pEnd*/, size_t numSymbols, size_t& numBitsDecoded, Functor& functor);

    // Encode a stream of characters
    inline size_t encodeStream(const CharDeque& input, const HuffmanTreeCodec<char>& symbolEncoder, HuffmanTreeCodec<int> & runEncoder, ByteVector& output)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(symbolEncoder.getMaxBits() <= BITS_PER_BYTE);
        assert(runEncoder.getMaxBits() <= BITS_PER_BYTE);

        // Re-construct the stream as a sequence of <symbol, run> pairs
        // Only allow runs whose length is encoded in the runEncoder
        CharDeque stream = input;
        RLEPairVector rlePairVector;
        while(!stream.empty())
        {
            // Count the current run
            char currChar = stream.front();
            size_t idx = 0;
            while(idx != stream.size() && stream[idx] == currChar)
                    ++idx;

            // A run of length idx has ended
            RLEPair pair(currChar, idx);
            rlePairVector.push_back(pair);

            // Extract idx characters from the stream
            for(size_t i = 0; i < idx; ++i)
                stream.pop_front();
        }

        // Count the number of bits it will require to encode the pairs
        size_t numBits = 0;
        for(RLEPairVector::const_iterator iter = rlePairVector.begin(); iter != rlePairVector.end(); ++iter)
        {
            numBits += 3;
            if(iter->second > 1)
                numBits += RLE::getEncodedLength(iter->second);
        }
        size_t numBytes = (numBits % BITS_PER_BYTE == 0) ? numBits / BITS_PER_BYTE : (numBits / BITS_PER_BYTE) + 1;
        
        //std::cout << numBits << " bits required to encode the stream\n";
        //std::cout << numBytes << " bytes required to encode the stream\n";

        output.resize(numBytes, 0);

        // Perform the encoding
        size_t currBit = 0;
        for(RLEPairVector::const_iterator iter = rlePairVector.begin(); iter != rlePairVector.end(); ++iter)
        {
#ifdef DEBUG_ENCODING
            std::cout << "Encoding pair: " << iter->first << "," << iter->second << "\n";
#endif
            currBit += _writeCode(BWT_ALPHABET::getRank(iter->first), 3, currBit, output);

            // Only encode the run if its length is greater than 1
            int rl = iter->second;
            if(rl > 1)
            {
                // Intialize the run length encoder
                RLE::initEncoding(rl);

                // Emit run codes until the entire run has been encoded
                int code = RLE::nextCode(rl);
                assert(code != 0);
                do
                {
                    currBit += _writeCode(code, RLE_SYMBOL_LEN, currBit, output);
                    code = RLE::nextCode(rl);
                } while(code != 0);
            }
        }

        std::string out;
        StringDecode sd(out);
        size_t nbd = 0;
        decodeStream(&(*output.begin()), &(*output.end()), input.size(), nbd, sd);

        std::string in(input.begin(), input.end());
//        std::cout << "IN:  " << in << "\n";
//        std::cout << "OUT: " << out << "\n";
        assert(in == out);
        return input.size();
    }

    // Decode a stream into the provided functor

#define DECODE_UNIT uint32_t
#define DECODE_UNIT_BYTES sizeof(DECODE_UNIT)
#define DECODE_UNIT_BITS DECODE_UNIT_BYTES * 8
#define DECODE_READ_LENGTH 3
#define DECODE_MASK 7
    // Decompress the data starting at pInput. The read cannot exceed the endpoint given by pEnd. Returns
    // the total number of symbols decoded. The out parameters numBitsDecoded is also set.
    template<typename Functor>
    inline size_t decodeStream(const unsigned char* pInput, const unsigned char* /*pEnd*/, size_t numSymbols, size_t& numBitsDecoded, Functor& functor)
    {
        size_t numSymbolsDecoded = 0;
        (void)numBitsDecoded;
        // Prime the decode unit by reading 16 bits from the stream
        DECODE_UNIT decodeUnit = pInput[0] << 24 | pInput[1] << 16 | pInput[2] << 8 | pInput[3];
        pInput += 4;
        DECODE_UNIT numBitsAvailable = 32;
        DECODE_UNIT code;

        // prime the read with the first symbol
        int currSym = (decodeUnit >> (numBitsAvailable - DECODE_READ_LENGTH)) & DECODE_MASK;
        numBitsAvailable -= DECODE_READ_LENGTH;
        int currRL = 0;
        int currRunIdx = 0;
        int target = numSymbols - numSymbolsDecoded;
        while(target > 0)
        {
            // Read in more bits if necessary
            if(numBitsAvailable < 8)
            {
                decodeUnit <<= 8;
                decodeUnit |= *pInput++;
                numBitsAvailable += 8;
            }

            code = (decodeUnit >> (numBitsAvailable - DECODE_READ_LENGTH) & DECODE_MASK);
            numBitsAvailable -= DECODE_READ_LENGTH;

            // Parse the code. If it was a symbol literal, add the current
            // run (if any) to the decoder. Otherwise, add the run symbol
            // to the current run. If the current run length matches
            // the target number of symbols to decode, stop the decoding
            if(code < RUN_A)
            {
                // If no run was encoded, the run length is 1
                currRL = currRunIdx == 0 ? 1 : currRL;
                functor(currSym, currRL);
                numSymbolsDecoded += currRL;
                target = numSymbols - numSymbolsDecoded;

                currSym = code;
                currRunIdx = 0;
                currRL = 0;
            }
            else
            {
                // Add the run code into the run length without branching
                // this calculates 2 ^ (currRunIdx) for RUN_A
                // and 2 ^ (currRunIdx + 1) for RUN_B
                currRL += (1 << (currRunIdx + (code - RUN_A)));
                ++currRunIdx;

                // If the target is reached, stop decoding
                if(currRL >= target)
                {
                    functor(currSym, target);
                    target = 0;
                }
            }
        }
        return numSymbols;
    }    
};

#endif
