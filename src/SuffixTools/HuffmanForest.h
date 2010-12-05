//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffmanForest - Singleton class containing a 
// family of huffman tree encoders. When encoding a set
// of symbols, the forest can be queried to return a tree
// that can efficiently encoded the data. When decoding,
// the appropriate decoder can be found by an index.
//
#ifndef HUFFMAN_FOREST_H
#define HUFFMAN_FOREST_H

#include "Huffman.h"

class HuffmanForest
{

    public:

        inline static HuffmanForest& Instance()
        {
            static HuffmanForest instance;
            return instance;
        }

        // Get an encoder for the data. outIdx is set to the index of 
        // the huffman tree used
        HuffmanTreeCodec<char>& getBestEncoder(const std::map<char, int>& symbolCounts, size_t& outIdx);

        inline HuffmanTreeCodec<char>& getDecoder(int idx)
        {
            return m_trees[idx];
        }

    private:
        HuffmanForest();
        ~HuffmanForest();

        //
        std::vector<HuffmanTreeCodec<char> > m_trees;

};

#endif
