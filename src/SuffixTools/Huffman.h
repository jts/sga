//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// Huffman - Simple huffman encoder/decoder for small alphabet
//
#ifndef HUFFMAN_H
#define HUFFMAN_H

#define MAX_BITS 4
#include <inttypes.h>
#include <map>
#include <vector>
#include <string>
#include "Alphabet.h"

struct HuffmanEncodePair
{
    uint8_t code;
    uint8_t bits;
};

struct HuffmanDecodePair
{
    char base;
    uint8_t bits;
};

typedef std::map<char, HuffmanEncodePair> HuffmanEncodeMap;
typedef std::vector<HuffmanDecodePair> HuffmanDecodeVector;

namespace Huffman
{
void buildHuffman(AlphaCount64& ac, HuffmanEncodeMap& outEncoder, HuffmanDecodeVector& outDecoder);
};
#endif
