//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BitVector - Vector of bits 
//
#include "BitVector.h"
#include <assert.h>

//
BitVector::BitVector()
{

}

//
BitVector::BitVector(size_t n)
{
    resize(n);
}

//
BitVector::~BitVector()
{

}

//
void BitVector::resize(size_t n)
{
    size_t num_bytes = (n % 8 == 0) ? n / 8 : (n / 8) + 1;
    m_data.resize(num_bytes);
}

// Set bit at position i to value v
void BitVector::set(size_t i, bool v)
{
    size_t byte = i / 8;
    assert(byte < m_data.size());
    size_t offset = i - byte * 8;
    m_data[byte].set(offset, v);
}

// Test bit i
bool BitVector::test(size_t i) const
{
    size_t byte = i / 8;
    size_t offset = i - byte * 8;
    return m_data[byte].test(offset);
}

