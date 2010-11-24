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

typedef std::map<char, HuffmanEncodePair> HuffmanSymbolEncoder;
typedef std::vector<HuffmanDecodePair> HuffmanSymbolDecoder;
typedef std::map<int, HuffmanEncodePair> HuffmanRunLengthEncoder;

// Huffman Tree
template<typename T>
struct HuffmanNode
{
    HuffmanNode(T sym, double f) : pParent(NULL), pLeftChild(NULL), pRightChild(NULL), symbol(sym), frequency(f) {}
    bool isLeaf() const { return pLeftChild == NULL && pRightChild == NULL; }
    HuffmanNode* pParent;
    HuffmanNode* pLeftChild;
    HuffmanNode* pRightChild;
    T symbol;
    double frequency;
};

// Returns true if pLHS has a greater frequency than pRHS
template<typename T>
struct HuffmanNodePriority
{
    bool operator()(HuffmanNode<T>* pLHS, HuffmanNode<T>* pRHS) const
    {
        return (pLHS->frequency > pRHS->frequency);
    }
};

namespace Huffman
{
void buildRLHuffmanTree(HuffmanRunLengthEncoder& outEncoder);
void visitHuffmanTree(HuffmanNode<int>* pNode, std::string code, HuffmanRunLengthEncoder& outEncoder);
int getEncodingLength(HuffmanRunLengthEncoder& outEncoder, int inLength);
void buildSymbolHuffman(AlphaCount64& ac, HuffmanSymbolEncoder& outEncoder, HuffmanSymbolDecoder& outDecoder);
void buildRunLengthHuffman(HuffmanRunLengthEncoder& outEncoder);

};

#endif
