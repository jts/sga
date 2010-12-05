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

#define BITS_MASK 255
#define PACKED_DECODE_SHIFT 8
#define PACKED_DECODE_TYPE int

class PackedDecodeTable
{
    public:
        PackedDecodeTable() {}

        void initialize(const HuffmanTreeCodec<int>& tree)
        {
            size_t max = tree.getMaxCode();
            m_readLen = tree.getMaxBits();
            m_decodeTable.reserve(max+1);
            for(size_t i = 0; i <= max; ++i)
            {
                m_decodeTable.push_back(pack(tree.decodeSymbol(i), tree.decodeBits(i)));
            }
        }

        inline int getCodeReadLength() const
        {
            return m_readLen;
        }

        inline PACKED_DECODE_TYPE pack(PACKED_DECODE_TYPE symbol, PACKED_DECODE_TYPE bits)
        {
            return (symbol << PACKED_DECODE_SHIFT) | bits;
        }

        inline void unpack(int code, PACKED_DECODE_TYPE& symOut, PACKED_DECODE_TYPE& bitsOut) const
        {
            PACKED_DECODE_TYPE in = m_decodeTable[code];
            bitsOut = in & BITS_MASK;
            symOut = in >> PACKED_DECODE_SHIFT;
        }

        std::vector<PACKED_DECODE_TYPE> m_decodeTable;
        int m_readLen;
};

namespace Huffman
{

HuffmanTreeCodec<int> buildRLHuffmanTree();
void buildSymbolHuffman(AlphaCount64& ac, HuffmanSymbolEncoder& outEncoder);
void buildRunLengthHuffman(HuffmanRunLengthEncoder& outEncoder);

};

#endif
