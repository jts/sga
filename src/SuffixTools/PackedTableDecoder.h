//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// PackedTableDecoder - Single-table decoders for huffman
// encoded data.
//
#ifndef PACKEDTABLEDECODER_H
#define PACKEDTABLEDECODER_H

#define BITS_MASK 255
#define PACKED_DECODE_SHIFT 8
#define PACKED_DECODE_TYPE int

class RLPackedTableDecoder
{
    public:
        RLPackedTableDecoder() {}

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

// Packed table decoder for characters
class CharPackedTableDecoder
{
    public:
        CharPackedTableDecoder() {}

        void initialize(const HuffmanTreeCodec<char>& tree)
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

        inline PACKED_DECODE_TYPE pack(int symbol, PACKED_DECODE_TYPE bits)
        {
            return (BWT_ALPHABET::getRank(symbol) << PACKED_DECODE_SHIFT) | bits;
        }

        // Unpack the symbol and coding bits.
        // This returns the symbols rank in the bwt alphabet, not the actual symbol
        inline void unpack(int code, int& rankOut, PACKED_DECODE_TYPE& bitsOut) const
        {
            PACKED_DECODE_TYPE in = m_decodeTable[code];
            bitsOut = in & BITS_MASK;
            rankOut = in >> PACKED_DECODE_SHIFT;
        }

        std::vector<PACKED_DECODE_TYPE> m_decodeTable;
        int m_readLen;
};

#endif
