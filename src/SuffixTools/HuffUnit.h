//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffUnit - A huffman-encoded unit of the
//
//
#ifndef HUFFUNIT_H
#define HUFFUNIT_H

#include <algorithm>

//
#define HUFF_BASE_MASK 0x07  //0000 0000 0000 0111
#define MAX_SYMBOLS 5
#define HUFF_MAXLEN 3

// A HuffmanUnit is a string of symbols encoded
// in 16 bits
struct HuffUnit
{
    HuffUnit() : data(0) {}

    // Add this data to the AlphaCount
    // Only add up to maxCount symbols. Returns the number
    // of symbols added
    inline size_t addAlphaCount(AlphaCount64& ac, size_t max) const
    {
        size_t count = std::min(max, (size_t)MAX_SYMBOLS);
        for(size_t i = 0; i < count; ++i)
            ac.increment(getChar(i));
        return count;
    }

    // Add this run to the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max
    inline size_t addCount(char b, size_t& base_count, size_t max) const
    {
        size_t count = std::min(max, (size_t)MAX_SYMBOLS);
        for(size_t i = 0; i < count; ++i)
        {
            char symbol = getChar(i);
            if(symbol == b)
                ++base_count;
        }
        return count;        
    }    

    // Subtract this unit from the AlphaCount
    // Only subtract up to maxCount symbols. Returns the number
    // of symbols subtracted.
    // Since this function is used to count backwards through a BWT,
    // we start from the back of the unit.
    inline size_t subtractAlphaCount(AlphaCount64& ac, size_t max) const
    {
        size_t count = std::min(max, (size_t)MAX_SYMBOLS);
        int stop = 4 - count;
        for(int i = 4; i > stop; --i)
            ac.decrement(getChar(i));
        return count;
    }

    // Subtract this run from the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max.
    // Since this function is used to count backwards through a BWT,
    // we start from the back of the unit.
    inline size_t subtractCount(char b, size_t& base_count, size_t max) const
    {
        size_t count = std::min(max, (size_t)MAX_SYMBOLS);
        int stop = 4 - count;
        for(int i = 4; i > stop; --i)
        {
            char symbol = getChar(i);
            if(symbol == b)
                --base_count;
        }
        return count;
    }
        
    // Set the symbol
    inline void setChar(char symbol, size_t idx)
    {
        // Shift mask
        uint16_t shift = (MAX_SYMBOLS - 1 - idx) * HUFF_MAXLEN;
        uint16_t mask = ~(HUFF_BASE_MASK << shift);
        
        // Clear current symbol
        data &= mask;
        
        uint16_t code = BWT_ALPHABET::getRank(symbol);
        code <<= shift;
        data |= code;
    }

    // Get the symbol
    inline char getChar(size_t idx) const
    {
        uint16_t shift = (MAX_SYMBOLS - 1 - idx) * HUFF_MAXLEN;
        uint16_t mask = (HUFF_BASE_MASK << shift);
        uint16_t code = data & mask;
        code >>= shift;
        return BWT_ALPHABET::getChar(code);
    }

    void print() const
    {
        for(int i = 0; i < MAX_SYMBOLS; ++i)
        {
            std::cout << getChar(i);
        }
        std::cout << "\n";
    }

    // 
    uint16_t data;

    friend class RLBWTReader;
    friend class RLBWTWriter;
};
typedef std::vector<RLUnit> RLVector;

#endif
