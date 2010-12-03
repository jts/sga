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

#include "Huffman.h"

//#define DEBUG_ENCODING 1
#define BITS_PER_BYTE 8

typedef std::pair<char, int> RLEPair;
typedef std::vector<RLEPair> RLEPairVector;
typedef std::deque<char> CharDeque;
typedef std::vector<uint8_t> EncodedArray;

namespace StreamEncode
{

    // Decode functors for the generic decoding function
    struct AlphaCountDecode
    {
        AlphaCountDecode(AlphaCount64& target) : m_target(target) {}
        inline void operator()(char b, int rl)
        {
            m_target.add(b, rl);
        }
        AlphaCount64& m_target;
    };

    struct StringDecode
    {
        StringDecode(std::string& target) : m_target(target) {}
        inline void operator()(char b, int rl)
        {
            m_target.append(rl, b);
        }
        std::string& m_target;
    };
    
    struct BaseCountDecode
    {
        BaseCountDecode(char targetBase, size_t& targetCount) : m_targetBase(targetBase), 
                                                                    m_targetCount(targetCount) {}
        inline void operator()(char b, int rl)
        {
            if(b == m_targetBase)
                m_targetCount += rl;
        }
        char m_targetBase;
        size_t& m_targetCount;
    };   

    //
    inline void printEncoding(const EncodedArray& output)
    {
        std::cout << "Encoding: ";
        for(size_t i = 0; i < output.size(); ++i)
        {
            std::cout << int2Binary(output[i], 8) << " ";
        }
        std::cout << "\n";
    }

    // Write the code into the output stream starting at currBit
    inline size_t _writeCode(EncodePair& ep, size_t currBit, EncodedArray& output)
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
    inline size_t _readCode(size_t currBit, size_t maxBits, const unsigned char* pInput, int& outCode)
    {

#ifdef DEBUG_ENCODING
        std::cout << "Reading " << maxBits << " from array starting at " << currBit << "\n";
#endif        
        outCode = 0;
        int bitsRemaining = maxBits;
        while(bitsRemaining > 0)
        {
            int byte = currBit / BITS_PER_BYTE;
            int bitIdx = MOD_POWER_2(currBit, BITS_PER_BYTE);
            int bitsToRead = std::min((BITS_PER_BYTE - bitIdx), bitsRemaining);

            //
            int mask = ((1 << bitsToRead) - 1);
            int currPos = (BITS_PER_BYTE - bitsToRead);

#ifdef DEBUG_ENCODING
            std::cout << "Mask " << int2Binary(mask, bitsToRead) << "\n";
            std::cout << "CurrPos: " << currPos << "\n";
            std::cout << "BitIdx: " << bitIdx << "\n";
#endif  
            // Shift the mask into place
            if(currPos < bitIdx)
                mask >>= (bitIdx - currPos);
            else if(currPos > bitIdx)
                mask <<= (currPos - bitIdx);

#ifdef DEBUG_ENCODING
            std::cout << "Shift mask " << int2Binary(mask, 8) << "\n";
#endif           
            // Read the value
            int tmpCode = (pInput[byte] & mask) >> (BITS_PER_BYTE - bitIdx - bitsToRead);

#ifdef DEBUG_ENCODING
            std::cout << "TmpC " << int2Binary(tmpCode, 8) << "\n";
#endif
            // Add the value into the outCode
            outCode <<= bitsToRead;
            outCode |= tmpCode;

            currBit += bitsToRead;
            bitsRemaining -= bitsToRead;
        }
#ifdef DEBUG_ENCODING
        std::cout << "out: " << int2Binary(outCode, maxBits) << "\n";
#endif        
        return maxBits;
    }
    
    // Encode a stream of characters
    inline size_t encodeStream(const CharDeque& input, const HuffmanTreeCodec<char>& symbolEncoder, HuffmanTreeCodec<int> & runEncoder, EncodedArray& output)
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
            size_t encodingLength = 1;//runEncoder.getGreatestLowerBound(idx);
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
    template<typename Functor>
    inline size_t decodeStream(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const unsigned char* pInput, size_t numSymbols, Functor& functor)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(symbolEncoder.getMaxBits() <= BITS_PER_BYTE);
        assert(runEncoder.getMaxBits() <= BITS_PER_BYTE);
        
        size_t symbolReadLen = symbolEncoder.getMaxBits();
        size_t runReadLen = runEncoder.getMaxBits();

        size_t numBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        while(numSymbolsDecoded < numSymbols)
        {
            // Read a symbol then a run
            int code = 0;
            _readCode(numBitsDecoded, symbolReadLen, pInput, code);

            // Parse the code
            HuffmanTreeCodec<char>::DecodePair sdp = symbolEncoder.decode(code);
            numBitsDecoded += sdp.bits;

            (void)runReadLen;
            _readCode(numBitsDecoded, runReadLen, pInput, code);
            HuffmanTreeCodec<int>::DecodePair rdp = runEncoder.decode(code);

            numBitsDecoded += rdp.bits;
            size_t addLength = std::min(rdp.symbol, (int)numSymbols - (int)numSymbolsDecoded);
            numSymbolsDecoded += addLength;

#ifdef DEBUG_VALIDATE
            std::cout << "Decoded pair: " << sdp.symbol << "," << rdp.symbol << "\n";
#endif
            functor(sdp.symbol, addLength);
        }
        return numSymbols;
    }
};

#endif
