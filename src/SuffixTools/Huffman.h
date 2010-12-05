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
#include "HuffmanTreeCodec.h"

typedef std::map<char, EncodePair> HuffmanSymbolEncoder;
typedef std::map<int, EncodePair> HuffmanRunLengthEncoder;


namespace Huffman
{

HuffmanTreeCodec<int> buildRLHuffmanTree();
void buildSymbolHuffman(AlphaCount64& ac, HuffmanSymbolEncoder& outEncoder);
void buildRunLengthHuffman(HuffmanRunLengthEncoder& outEncoder);
};

#endif
