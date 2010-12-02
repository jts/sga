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
#define BITS_AVAILABLE 64

typedef std::pair<char, int> RLEPair;
typedef std::vector<RLEPair> RLEPairVector;
typedef std::deque<char> CharDeque;
typedef std::vector<uint8_t> EncodedArray;
typedef uint64_t EncodedValue;
typedef std::vector<EncodedValue> EVVector;

namespace StreamEncode
{
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

    // Decode an array into a string
    inline size_t decodeStream(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const unsigned char* pInput, size_t numSymbols, std::string& decoded)
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
            const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            numBitsDecoded += sdp.bits;

            _readCode(numBitsDecoded, runReadLen, pInput, code);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(code);
            numBitsDecoded += rdp.bits;
            numSymbolsDecoded += rdp.symbol;

#ifdef DEBUG_VALIDATE
            std::cout << "Decoded pair: " << sdp.symbol << "," << rdp.symbol << "\n";
#endif
            decoded.append(rdp.symbol, sdp.symbol);
        }
        return numSymbols;
    }

    // Decode from the array pInput up to numSymbols. Add the counts of symbols to counts
    inline size_t countDecoded(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const unsigned char* pInput, size_t numSymbols, AlphaCount64& counts)
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
            const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            numBitsDecoded += sdp.bits;

            _readCode(numBitsDecoded, runReadLen, pInput, code);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(code);
            numBitsDecoded += rdp.bits;
            
            // Cap the run length to add at the number of symbols left to process
            size_t addLength = std::min(rdp.symbol, (int)numSymbols - (int)numSymbolsDecoded);
            numSymbolsDecoded += addLength;
            counts.add(sdp.symbol, addLength);
        }
        return numSymbolsDecoded;
    }

        // Decode from the array pInput up to numSymbols. Add the counts of symbols to counts
    inline size_t countDecoded(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const unsigned char* pInput, size_t numSymbols, char b, size_t& count)
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
            const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            numBitsDecoded += sdp.bits;

            _readCode(numBitsDecoded, runReadLen, pInput, code);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(code);
            numBitsDecoded += rdp.bits;
            
            // Cap the run length to add at the number of symbols left to process
            size_t addLength = std::min(rdp.symbol, (int)numSymbols - (int)numSymbolsDecoded);
            numSymbolsDecoded += addLength;

            if(sdp.symbol == b)
                count += addLength;
        }
        return numSymbolsDecoded;
    }
};

namespace StreamEncode2
{
    //
    inline void printEncoding(const EncodedValue output)
    {
        std::cout << "Encoding: " << int2Binary(output, BITS_AVAILABLE) << "\n";
    }

    // Write the code into the output stream starting at currBit. Returns the number of bits written
    inline size_t _writeCode(EncodePair& ep, size_t currBit, EncodedValue& output)
    {
#ifdef DEBUG_ENCODING
        printEncoding(output);
        std::cout << "Writing the code " << int2Binary(ep.code, ep.bits) << " at bit " << currBit << "\n";
#endif
        size_t code = ep.code;
        int codeBits = ep.bits;

        // Calculate the shift values and masks to apply
        int shiftVal = BITS_AVAILABLE - currBit - codeBits;
        // Set the mask and shift it into position
        size_t mask = (1 << codeBits) - 1;
        size_t inCode = (code & mask);
        
#ifdef DEBUG_ENCODING
        std::cout << "Mask: " << int2Binary(mask, codeBits) << "\n";
        std::cout << "Masked: " << int2Binary(inCode,codeBits) << "\n";
#endif
        inCode <<= shiftVal;
#ifdef DEBUG_ENCODING
        std::cout << "inCode: " << int2Binary(inCode, BITS_AVAILABLE) << "\n";
#endif 
        output |= inCode;
        return codeBits;
    }

    // Read maxBits from the array starting at currBits and write the value to outCode
    // Returns the number of bits read
    inline size_t _readCode(size_t currBit, size_t maxBits, const EncodedValue& input, EncodedValue& outCode)
    {

#ifdef DEBUG_ENCODING
        printEncoding(input);
        std::cout << "Reading " << maxBits << " from array starting at " << currBit << "\n";
#endif
        outCode = 0;
        // If maxBits is larger than the number of bits left in the unit, we
        // only read bitsRemaining and shift the code up by the difference. This is
        // to handle the case where a run was encoded at the end of a unit that is less than 8
        // bits in length. The bits that were not directly read are implicitly assumed to be zero
        size_t bitsRemaining = BITS_AVAILABLE - currBit;

        if(bitsRemaining > maxBits)
        {
            size_t mask = ((1 << maxBits) - 1);
            int shiftVal = (bitsRemaining - maxBits);

#ifdef DEBUG_ENCODING
            std::cout << "Mask " << int2Binary(mask, maxBits) << "\n";
#endif      
                   
            mask <<= shiftVal;

#ifdef DEBUG_ENCODING
            std::cout << "ShiftMask " << int2Binary(mask, maxBits) << "\n";
#endif 

            outCode = (input & mask) >> shiftVal;

#ifdef DEBUG_ENCODING
            std::cout << "Outcode: " << int2Binary(outCode, maxBits) << "\n";
#endif
            return maxBits;
        }
        else
        {
            size_t mask = ((1 << bitsRemaining) - 1);

#ifdef DEBUG_ENCODING
            std::cout << "Mask " << int2Binary(mask, bitsRemaining) << "\n";
#endif                
            int diff = maxBits - bitsRemaining;

            outCode = (input & mask) << diff;
            return maxBits;
        }
        return maxBits;
    }
    
    // Encode a stream of characters
    inline size_t encodeStream(const CharDeque& input, const HuffmanTreeCodec<char>& symbolEncoder, HuffmanTreeCodec<int> & runEncoder, EVVector& outputArray)
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

        outputArray.clear();
        // Perform the encoding
        size_t currBit = 0;
        EncodedValue currentValue = 0;
        EncodePair terminationCode = symbolEncoder.encode('\0');

        RLEPairVector::const_iterator iter = rlePairVector.begin();
        while(iter != rlePairVector.end())
        {
#ifdef DEBUG_ENCODING
            std::cout << "Encoding pair: " << iter->first << "," << iter->second << "\n";
#endif
            EncodePair symEP = symbolEncoder.encode(iter->first);
            EncodePair rlEP = runEncoder.encode(iter->second);

            if(BITS_AVAILABLE - currBit - terminationCode.bits > symEP.bits + rlEP.bits)
            {
                currBit += _writeCode(symEP, currBit, currentValue);
                currBit += _writeCode(rlEP, currBit, currentValue);
                ++iter;
              
            }
            else
            {
                // Add a termination symbol to the value
                //std::cout << "writing termination code!\n";
                assert(BITS_AVAILABLE - currBit > terminationCode.bits);
                currBit += _writeCode(terminationCode, currBit, currentValue);
                outputArray.push_back(currentValue);
                currentValue = 0;
                currBit = 0;
            }
        }
        
        // Terminate the last code if it has any data and push it
        if(currBit > 0)
        {
            assert(BITS_AVAILABLE - currBit > terminationCode.bits);
            currBit += _writeCode(terminationCode, currBit, currentValue);
            outputArray.push_back(currentValue);
        }
        return input.size();
    }

    // Decode an encodedarray into a string
    inline size_t decodeStream(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const EncodedValue* pInput, size_t numSymbols, std::string& decoded)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(symbolEncoder.getMaxBits() <= BITS_PER_BYTE);
        assert(runEncoder.getMaxBits() <= BITS_PER_BYTE);
        
        size_t symbolReadLen = symbolEncoder.getMaxBits();
        size_t runReadLen = runEncoder.getMaxBits();

        size_t currBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        while(numSymbolsDecoded < numSymbols)
        {
            // Read a symbol then a run
            EncodedValue code = 0;
            _readCode(currBitsDecoded, symbolReadLen, *pInput, code);

            // Parse the code
            const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            currBitsDecoded += sdp.bits;

            if(sdp.symbol == '\0')
            {
                // termination code, go to the next unit
                pInput += 1;
                currBitsDecoded = 0;
                continue;
            }   

            _readCode(currBitsDecoded, runReadLen, *pInput, code);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(code);

            currBitsDecoded += rdp.bits;
            numSymbolsDecoded += rdp.symbol;

#ifdef DEBUG_ENCODING
            std::cout << "Decoded pair: " << sdp.symbol << "," << rdp.symbol << "\n";
#endif
            decoded.append(rdp.symbol, sdp.symbol);
#ifdef DEBUG_ENCODING
            std::cout << "Decoding so far: " << decoded << "\n";
#endif
        }
        return numSymbols;
    }

    // Decode from the array pInput up to numSymbols. Add the counts of symbols to counts
    inline size_t countDecoded(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const EncodedValue* pInput, size_t numSymbols, AlphaCount64& counts)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(symbolEncoder.getMaxBits() <= BITS_PER_BYTE);
        assert(runEncoder.getMaxBits() <= BITS_PER_BYTE);
        
        size_t symbolReadLen = symbolEncoder.getMaxBits();
        size_t runReadLen = runEncoder.getMaxBits();

        size_t currBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        while(numSymbolsDecoded < numSymbols)
        {
            // Read a symbol then a run
            EncodedValue symCode = 0;
            _readCode(currBitsDecoded, symbolReadLen, *pInput, symCode);

            // Parse the code
            //const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(symCode);
            currBitsDecoded += 3;//sdp.bits;

            if(symCode == 7)
            {
                // termination code, go to the next unit
                pInput += 1;
                currBitsDecoded = 0;
                continue;
            }

            EncodedValue rlCode = 0;
            _readCode(currBitsDecoded, runReadLen, *pInput, rlCode);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(rlCode);
            int rl = rdp.symbol;
            currBitsDecoded += rdp.bits;
            
            // Cap the run length to add at the number of symbols left to process
            size_t addLength = std::min(rl, (int)numSymbols - (int)numSymbolsDecoded);
            numSymbolsDecoded += addLength;
            counts.add(ALPHABET[symCode], addLength);
        }
        return numSymbolsDecoded;
    }

    // Decode from the array pInput up to numSymbols. Add the counts of symbols to counts
    inline size_t countDecoded(const HuffmanTreeCodec<char>& symbolEncoder, const HuffmanTreeCodec<int>& runEncoder, const EncodedValue* pInput, size_t numSymbols, char b, size_t& count)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(symbolEncoder.getMaxBits() <= BITS_PER_BYTE);
        assert(runEncoder.getMaxBits() <= BITS_PER_BYTE);
        
        size_t symbolReadLen = symbolEncoder.getMaxBits();
        size_t runReadLen = runEncoder.getMaxBits();

        size_t currBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        while(numSymbolsDecoded < numSymbols)
        {
            // Read a symbol then a run
            EncodedValue code = 0;
            _readCode(currBitsDecoded, symbolReadLen, *pInput, code);

            // Parse the code
            const HuffmanTreeCodec<char>::DecodePair& sdp = symbolEncoder.decode(code);
            currBitsDecoded += sdp.bits;

            if(sdp.symbol == '\0')
            {
                // termination code, go to the next unit
                pInput += 1;
                currBitsDecoded = 0;
                continue;
            }

            _readCode(currBitsDecoded, runReadLen, *pInput, code);
            const HuffmanTreeCodec<int>::DecodePair& rdp = runEncoder.decode(code);
            currBitsDecoded += rdp.bits;
            
            // Cap the run length to add at the number of symbols left to process
            size_t addLength = std::min(rdp.symbol, (int)numSymbols - (int)numSymbolsDecoded);
            numSymbolsDecoded += addLength;

            if(sdp.symbol == b)
                count += addLength;
        }
        return numSymbolsDecoded;
    }
};

#endif
