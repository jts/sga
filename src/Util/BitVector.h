//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BitVector - Vector of bits 
//
#ifndef BITVECTOR_H
#define BITVECTOR_H

#include "BitChar.h"
#include <vector>

class BitVector
{
    public:
    
        BitVector();
        BitVector(size_t n);
        ~BitVector();

        void resize(size_t n);
        void set(size_t i, bool v);
        bool test(size_t i) const;

    private:
        std::vector<BitChar> m_data;

};

#endif
