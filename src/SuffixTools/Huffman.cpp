//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// Huffman - Simple huffman encoder/decoder for small alphabet
//
#include "Huffman.h"
#include "Alphabet.h"
#include <assert.h>
#include <iostream>

void Huffman::buildHuffman(AlphaCount64& counts, HuffmanEncodeMap& outEncoder, HuffmanDecodeVector& outDecoder)
{
    std::string sortString = counts.getSortString();

    // Count number of non-zero symbols
    size_t nonzero = 0;
    for(size_t i = 0; i < sortString.size(); ++i)
    {
        if(counts.get(sortString[i]) > 0)
            nonzero += 1;
    }

    uint8_t standard5[] = {0, 2, 6, 14, 15};
    uint8_t bits5[] = {1, 2, 3, 4, 4};
    //uint8_t standard5[] = {0, 4, 5, 6, 14};
    //uint8_t bits5[] = {1, 3, 3, 3, 4};

    uint8_t standard4[] = {0, 1, 2, 3, 15};
    uint8_t bits4[] = {2, 2, 2, 2, 4};

    uint8_t standard3[] = {0, 1, 3, 15, 15};
    uint8_t bits3[] = {1, 2, 2, 4, 4};
    
    uint8_t standard2[] = {0, 1, 15, 15, 15};
    uint8_t bits2[] = {1, 1, 4, 4, 4};

    uint8_t* pCodes;
    uint8_t* pBits;
    std::string codeStr;
    if(nonzero == 5)
    {
        pCodes = standard5;
        pBits = bits5;
    }
    else if(nonzero == 4)
    {
        pCodes = standard4;
        pBits = bits4;
    }
    else if(nonzero == 3)
    {
        pCodes = standard3;
        pBits = bits3;
    }
    else
    {
        pCodes = standard2;
        pBits = bits2;
    }

    size_t numBitsNeeded = 0;
    std::string blockStr;
    assert(sortString.size() <= 5);
    for(size_t i = 0; i < sortString.size(); ++i)
    {
        char b = sortString[i];
        HuffmanEncodePair ep = {pCodes[i], pBits[i]};
        outEncoder.insert(std::make_pair(b, ep));

        numBitsNeeded += pBits[i] * counts.get(b);
        blockStr.append(counts.get(b), b);
    }

    //size_t n = blockStr.size();
    //std::cout << "Bits needed for block of len " << blockStr.size() << " " << counts << " = " << numBitsNeeded << " " << (double)numBitsNeeded / n << "\n";
    (void)outDecoder;
}

