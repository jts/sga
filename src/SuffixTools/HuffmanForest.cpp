//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffmanForst - Set of huffman trees
//
#include "HuffmanForest.h"
#include <set>

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
    
    // 4-symbol code alphabet
    // 0, 10, 110, 111, N/A
    CodeVector fourSymCode;
    fourSymCode.push_back(std::make_pair(0, 1));
    fourSymCode.push_back(std::make_pair(2, 2));
    fourSymCode.push_back(std::make_pair(6, 3));
    fourSymCode.push_back(std::make_pair(7, 3));

    // 3-symbol code alphabet 
    // 0, 10, 11
    CodeVector threeSymCode;
    threeSymCode.push_back(std::make_pair(0, 1));
    threeSymCode.push_back(std::make_pair(2, 2));
    threeSymCode.push_back(std::make_pair(3, 2));

    std::vector<CodeVector*> codeVectorPtrVector;
    codeVectorPtrVector.push_back(&fiveSymCode);
    codeVectorPtrVector.push_back(&fourSymCode);
    //codeVectorPtrVector.push_back(&threeSymCode);

    /*
    fiveSymCode.push_back(std::make_pair(0, 3));
    fiveSymCode.push_back(std::make_pair(1, 3));
    fiveSymCode.push_back(std::make_pair(2, 3));
    fiveSymCode.push_back(std::make_pair(3, 3));
    fiveSymCode.push_back(std::make_pair(4, 3));
    */

    // For every permutation of ACGT$, construct a huffman tree for all codings
    for(std::vector<CodeVector*>::iterator iter = codeVectorPtrVector.begin();
        iter != codeVectorPtrVector.end(); ++iter)
    {
        const CodeVector& codeVector = *(*iter);

        // Calculate all the permutations of the string with a prefix of length codeVector.size()
        std::string symbols = "$ACGT";
        std::set<std::string> usedPermutations;
        while(std::next_permutation(symbols.begin(), symbols.end()))
        {
            std::string currPermutation = symbols.substr(0, codeVector.size());
            if(usedPermutations.find(currPermutation) == usedPermutations.end())
            {
                HuffmanTreeCodec<char> currTree;
                
                // Explicitly set the huffman tree codes
                for(size_t i = 0; i < codeVector.size(); ++i)
                    currTree.explicitCode(currPermutation[i], codeVector[i].first, codeVector[i].second);
                m_trees.push_back(currTree);

                // Construct decoder
                CharPackedTableDecoder decoder;
                decoder.initialize(currTree);
                m_decoders.push_back(decoder);

                usedPermutations.insert(currPermutation);
            }
        }
    }

    std::cout << "Huffman family has " << m_decoders.size() << " trees\n";
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
    assert(minBits != std::numeric_limits<size_t>::max());
    outIdx = minIdx;
    return m_trees[outIdx];
}
