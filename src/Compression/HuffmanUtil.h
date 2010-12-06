//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffmanUtil - Utility functions for huffman encoding
//
#ifndef HUFFMANUTIL_H
#define HUFFMANUTIL_H

#include <inttypes.h>
#include <map>
#include <vector>
#include <string>
#include "Alphabet.h"
#include "HuffmanTreeCodec.h"

namespace HuffmanUtil
{

HuffmanTreeCodec<int> buildRLHuffmanTree();

};

#endif
