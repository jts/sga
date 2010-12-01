//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffmanForst - Set of huffman trees
//
#include "HuffmanForest.h"

//
HuffmanForest::HuffmanForest()
{
    // Construct the default trees
    std::map<char, int> countMap;
    countMap['A'] = 1;
    countMap['C'] = 1;
    countMap['G'] = 1;
    countMap['T'] = 1;
    countMap['$'] = 1;

    HuffmanTreeCodec<char> defaultTree(countMap);
    m_trees.push_back(defaultTree);
}

//
HuffmanForest::~HuffmanForest()
{

}

//
HuffmanTreeCodec<char>& HuffmanForest::getEncoder(size_t& outIdx)
{
    outIdx = 0;
    return m_trees.back();
}
