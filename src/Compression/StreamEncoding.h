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
    inline size_t _writeCode(EncodePair& ep, size_t currBit, ByteVector& output)
    {
#ifdef DEBUG_ENCODING
        printEncoding(output);
        std::cout << "Writing the code " << int2Binary(ep.code, ep.bits) << " at bit " << currBit << "\n";
#endif

        size_t code = ep.code;
        int codeBits = ep.bits;
        int bitsRemaining = codeBits;
        int bitOffset = 0;

        while(bitsRemaining > 0)
        {
            // Calculate position to start the write
            int byte = currBit / BITS_PER_BYTE;
            int bitIdx = MOD_POWER_2(currBit, BITS_PER_BYTE);
            int bitsToWrite = std::min((BITS_PER_BYTE - bitIdx), bitsRemaining);

            // Calculate the shift values and masks to apply
            int currPos = (BITS_PER_BYTE - codeBits + bitOffset);
            
            // Mask off the bits we want
            int mask = ((1 << bitsToWrite) - 1) << (codeBits - (bitsToWrite + bitOffset));
            int inCode = code & mask;

#ifdef DEBUG_ENCODING
            std::cout << "Mask: " << int2Binary(mask, codeBits) << "\n";
            std::cout << "Masked: " << int2Binary(inCode,codeBits) << "\n";
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
        return codeBits;
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
            size_t encodingLength = runEncoder.getGreatestLowerBound(idx);
            RLEPair pair(currChar, encodingLength);
            rlePairVector.push_back(pair);

            // Extract idx characters from the stream
            for(size_t i = 0; i < encodingLength; ++i)
                stream.pop_front();
        }

        // Count the number of bits it will require to encode the pairs
        size_t numBits = 0;
        for(RLEPairVector::const_iterator iter = rlePairVector.begin(); iter != rlePairVector.end(); ++iter)
        {
            numBits += symbolEncoder.encode(iter->first).bits;
            numBits += runEncoder.encode(iter->second).bits;
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
            EncodePair symEP = symbolEncoder.encode(iter->first);
            currBit += _writeCode(symEP, currBit, output);
            EncodePair rlEP = runEncoder.encode(iter->second);
            currBit += _writeCode(rlEP, currBit, output);
        }

        return input.size();
    }

    // Decode a stream into the provided functor

#define DECODE_UNIT uint64_t
#define DECODE_UNIT_BYTES sizeof(DECODE_UNIT)
#define DECODE_UNIT_BITS DECODE_UNIT_BYTES * 8
    template<typename Functor>
    inline size_t decodeStream(const CharPackedTableDecoder& symbolDecoder, const RLPackedTableDecoder& rlDecoder, 
                               const unsigned char* pInput, size_t numSymbols, Functor& functor)
    {
        DECODE_UNIT symbolReadLen = symbolDecoder.getCodeReadLength();
        DECODE_UNIT runReadLen = rlDecoder.getCodeReadLength();

        DECODE_UNIT numBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        
        // Prime the decode unit by reading 16 bits from the stream
        DECODE_UNIT decodeUnit = *pInput++;
        for(size_t i = 0; i < DECODE_UNIT_BYTES - 1; ++i)
        {
            decodeUnit <<= BITS_PER_BYTE;
            decodeUnit |= *pInput++;
        }

        DECODE_UNIT symMask = (1 << symbolReadLen) - 1;
        DECODE_UNIT symBaseShift = DECODE_UNIT_BITS - symbolReadLen;
        DECODE_UNIT rlMask = (1 << runReadLen) - 1;
        DECODE_UNIT rlBaseShift = DECODE_UNIT_BITS - runReadLen;

        const std::vector<PACKED_DECODE_TYPE>* pDecodeChar = symbolDecoder.getTable();
        const std::vector<PACKED_DECODE_TYPE>* pDecodeRL = rlDecoder.getTable();

        while(1)
        {
            // Read a symbol then a run
            DECODE_UNIT code = 0;
            code = (decodeUnit >> (symBaseShift - numBitsDecoded)) & symMask;
            // Parse the code
            //const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            PACKED_DECODE_TYPE packedCode = (*pDecodeChar)[code];
            int symRank = UNPACK_SYMBOL(packedCode);
            numBitsDecoded += UNPACK_BITS(packedCode);

            code = (decodeUnit >> (rlBaseShift - numBitsDecoded)) & rlMask;
            packedCode = (*pDecodeRL)[code];
            int rl = UNPACK_SYMBOL(packedCode);
            numBitsDecoded += UNPACK_BITS(packedCode);

            int diff = numSymbols - numSymbolsDecoded;
            if(rl < diff)
            {
                functor(symRank, rl);
                numSymbolsDecoded += rl;
            }
            else
            {
                functor(symRank, diff);
                return numSymbolsDecoded + diff;
            }

            // Update the decode unit
            if(DECODE_UNIT_BITS - numBitsDecoded < 2*BITS_PER_BYTE)
            {
                for(size_t i = 2; i < DECODE_UNIT_BYTES; ++i)
                    decodeUnit = (decodeUnit << BITS_PER_BYTE) | *pInput++;
                numBitsDecoded -= (BITS_PER_BYTE * (DECODE_UNIT_BYTES - 2));
            }

        }
        return numSymbols;
    }    
};

#endif
