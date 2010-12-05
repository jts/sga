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
    // Construct the possible huffman codes over a 5-symbol alphabet
    // 0, 10, 110, 1110, 1111

    typedef std::vector<std::pair<int, int> > CodeVector;
    CodeVector fiveSymCode;
    
    fiveSymCode.push_back(std::make_pair(0, 1));
    fiveSymCode.push_back(std::make_pair(2, 2));
    fiveSymCode.push_back(std::make_pair(6, 3));
    fiveSymCode.push_back(std::make_pair(14, 4));
    fiveSymCode.push_back(std::make_pair(15, 4));
    
    /*
    fiveSymCode.push_back(std::make_pair(0, 3));
    fiveSymCode.push_back(std::make_pair(1, 3));
    fiveSymCode.push_back(std::make_pair(2, 3));
    fiveSymCode.push_back(std::make_pair(3, 3));
    fiveSymCode.push_back(std::make_pair(4, 3));
    */

    // For every permutation of ACGT$, construct a huffman tree
    std::string symbols = "$ACGT";

    while(std::next_permutation(symbols.begin(), symbols.end()))
    {
        HuffmanTreeCodec<char> currTree;
        
        // Explicitly set the huffman tree codes
        for(size_t i = 0; i < symbols.size(); ++i)
            currTree.explicitCode(symbols[i], fiveSymCode[i].first, fiveSymCode[i].second);
        m_trees.push_back(currTree);
    }
}

//
HuffmanForest::~HuffmanForest()
{

}

//
HuffmanTreeCodec<char>& HuffmanForest::getBestEncoder(const std::map<char, int>& symbolCounts, size_t& outIdx)
{
    size_t minBits = std::numeric_limits<size_t>::max();
    size_t minIdx = 0;
    for(size_t i = 0; i < m_trees.size(); ++i)
    {
        size_t requiredBits = m_trees[i].getRequiredBits(symbolCounts);
        if(requiredBits < minBits)
        {
            minBits = requiredBits;
            minIdx = i;
        }
    }
    outIdx = minIdx;
    return m_trees[outIdx];
}
