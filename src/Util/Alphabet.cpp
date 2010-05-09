//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Alphabet.h - Abstraction of the alphabet
// that is used in the suffix array, fm-index,
// etc data structures
//
#include "Alphabet.h"
#include <iostream>
#include <iterator>

//
std::ostream& operator<<(std::ostream& out, const AlphaCount& ac)
{
    std::copy(ac.m_counts, ac.m_counts+ALPHABET_SIZE, std::ostream_iterator<BaseCount>(out, " "));
    return out;
}

std::istream& operator>>(std::istream& in, AlphaCount& ac)
{
    for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        in >> ac.m_counts[i];
    return in;
}

