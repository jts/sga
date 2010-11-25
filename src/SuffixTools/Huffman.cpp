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
#include <queue>

void Huffman::buildSymbolHuffman(AlphaCount64& counts, HuffmanSymbolEncoder& outEncoder)
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

    //uint8_t standard4[] = {0, 2, 6, 7, 15};
    //uint8_t bits4[] = {1, 2, 3, 3, 4};

    uint8_t standard3[] = {0, 1, 3, 15, 15};
    uint8_t bits3[] = {1, 2, 2, 4, 4};
    
    uint8_t standard2[] = {0, 1, 15, 15, 15};
    uint8_t bits2[] = {1, 1, 4, 4, 4};

    uint8_t* pCodes;
    uint8_t* pBits;
    if(nonzero == 5)
    {
        pCodes = standard5;
        pBits = bits5;
        minBits = 4;
    }
    else if(nonzero == 4)
    {
        pCodes = standard4;
        pBits = bits4;
        minBits = 2;
    }
    else if(nonzero == 3)
    {
        pCodes = standard3;
        pBits = bits3;
        minBits = 2;
    }
    else
    {
        pCodes = standard2;
        pBits = bits2;
        minBits = 1;
    }
        
    //pCodes = standard5;
    //pBits = bits5;

    size_t numBitsNeeded = 0;
    std::string blockStr;
    assert(sortString.size() <= 5);
    for(size_t i = 0; i < sortString.size(); ++i)
    {
        char b = sortString[i];
        EncodePair ep = {pCodes[i], pBits[i]};
        outEncoder.insert(std::make_pair(b, ep));

        numBitsNeeded += pBits[i] * counts.get(b);
        blockStr.append(counts.get(b), b);
    }

    //size_t n = blockStr.size();
    //std::cout << "Bits needed for block of len " << blockStr.size() << " " << counts << " = " << numBitsNeeded << " " << (double)numBitsNeeded / n << "\n";
}

void Huffman::buildRunLengthHuffman(HuffmanRunLengthEncoder& outEncoder)
{
    //uint8_t rl[] =  {1, 2, 4,  8,  16, 32, 64};
    //uint8_t standard[] = {0, 2, 6, 14, 30, -1, -1};
    //uint8_t bits[] =     {1, 2, 3, 4, 5, 6, 6};
    
    uint8_t rl[] = {1, 2, 4, 8, 16, 32, 64};
    uint8_t standard[] = {0, 2, 6, 14, 30, -1, -1};
    uint8_t bits[] =     {1, 2, 3, 4, 5, 5, 6};

    for(size_t i = 0; i < 7; ++i)
    {
        EncodePair ep = {standard[i], bits[i]};
        outEncoder.insert(std::make_pair(rl[i], ep));
    }
}

HuffmanTreeCodec<int> Huffman::buildRLHuffmanTree()
{
    // Hardcoded counts 
    typedef std::map<int, int> IntIntMap;
    IntIntMap input;

    input.insert(std::make_pair(1,32663330));
    input.insert(std::make_pair(2,6627699));
    input.insert(std::make_pair(3,1889314));
    input.insert(std::make_pair(4,836362));
    input.insert(std::make_pair(5,746639));
    input.insert(std::make_pair(6,991933));
    input.insert(std::make_pair(7,1353892));
    input.insert(std::make_pair(8,1705611));
    input.insert(std::make_pair(9,1942392));
    input.insert(std::make_pair(10,2001377));
    input.insert(std::make_pair(11,1887666));
    input.insert(std::make_pair(12,1646809));
    input.insert(std::make_pair(13,1332408));
    input.insert(std::make_pair(14,1022598));
    input.insert(std::make_pair(15,756765));
    input.insert(std::make_pair(16,559331));
    /*input.insert(std::make_pair(17,428257));
    input.insert(std::make_pair(18,345921));
    input.insert(std::make_pair(19,297472));
    input.insert(std::make_pair(20,268562));
    input.insert(std::make_pair(21,247393));
    input.insert(std::make_pair(22,226436));
    input.insert(std::make_pair(23,204497));
    input.insert(std::make_pair(24,180382));
    input.insert(std::make_pair(25,156150));
    input.insert(std::make_pair(26,133070));
    input.insert(std::make_pair(27,112234));
    input.insert(std::make_pair(28,94537));
    input.insert(std::make_pair(29,80767));
    input.insert(std::make_pair(30,69259));
    input.insert(std::make_pair(31,60695));
    */
    input.insert(std::make_pair(32,52863));
    /*
    input.insert(std::make_pair(33,46801));
    input.insert(std::make_pair(34,42216));
    input.insert(std::make_pair(35,37341));
    input.insert(std::make_pair(36,33860));
    input.insert(std::make_pair(37,29999));
    input.insert(std::make_pair(38,27147));
    input.insert(std::make_pair(39,24767));
    input.insert(std::make_pair(40,22479));
    input.insert(std::make_pair(41,20439));
    input.insert(std::make_pair(42,18672));
    input.insert(std::make_pair(43,17603));
    input.insert(std::make_pair(44,16201));
    input.insert(std::make_pair(45,15388));
    input.insert(std::make_pair(46,14724));
    input.insert(std::make_pair(47,13855));
    input.insert(std::make_pair(48,13664));
    input.insert(std::make_pair(49,12583));
    input.insert(std::make_pair(50,11960));
    input.insert(std::make_pair(51,11216));
    input.insert(std::make_pair(52,10682));
    input.insert(std::make_pair(53,10043));
    input.insert(std::make_pair(54,9077));
    input.insert(std::make_pair(55,8717));
    input.insert(std::make_pair(56,8254));
    input.insert(std::make_pair(57,7571));
    input.insert(std::make_pair(58,6973));
    input.insert(std::make_pair(59,6473));
    input.insert(std::make_pair(60,6000));
    input.insert(std::make_pair(61,5529));
    input.insert(std::make_pair(62,4989));
    input.insert(std::make_pair(63,4808));
    */
    input.insert(std::make_pair(64,78833));
    return HuffmanTreeCodec<int>(input);
}
