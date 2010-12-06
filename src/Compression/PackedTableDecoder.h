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

#define PACKED_BITS_MASK 255
#define PACKED_DECODE_SHIFT 8
#define PACKED_DECODE_TYPE int

#define UNPACK_SYMBOL(in) (in) >> PACKED_DECODE_SHIFT
#define UNPACK_BITS(in) (in) & PACKED_BITS_MASK

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

        // Pack a value into the table
        inline PACKED_DECODE_TYPE pack(PACKED_DECODE_TYPE symbol, PACKED_DECODE_TYPE bits)
        {
            return (symbol << PACKED_DECODE_SHIFT) | bits;
        }

        // Return a pointer to the table
        inline const std::vector<PACKED_DECODE_TYPE>* getTable() const
        {
            return &m_decodeTable;
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

        //
        inline PACKED_DECODE_TYPE pack(int symbol, PACKED_DECODE_TYPE bits)
        {
            return (BWT_ALPHABET::getRank(symbol) << PACKED_DECODE_SHIFT) | bits;
        }

        // Return a pointer to the table
        inline const std::vector<PACKED_DECODE_TYPE>* getTable() const
        {
            return &m_decodeTable;
        }

        std::vector<PACKED_DECODE_TYPE> m_decodeTable;
        int m_readLen;
};

#endif
