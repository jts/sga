//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLE -- Run length encoding. Encode a non-zero 
// integer x as a stream of 3 bit symbols corresponding to 
// the binary representation of (x - 1). This is the 
// way bzip2 encodes runs of symbols
//
#ifndef RLE_H
#define RLB_H

// The symbols encoding the run
#define RUN_A 6
#define RUN_B 7
#define RLE_SYMBOL_LEN 3

namespace RLE
{

// Returns the number of bits required to encode
// a run of length x.
inline size_t getEncodedLength(int x)
{
    assert(x > 0);
    x += 1;
    size_t numSymbols = 0;
    while(x >>= 1)
    {
        numSymbols += 1;
    }
    return RLE_SYMBOL_LEN * numSymbols;
}

// Initialize the RLE encoder for a run of length x
inline void initEncoding(int& x)
{
    x -= 1;
}

// Get the next run code that should be output
// to the stream. If 0 is returned, stop the encoding of
// this run. x is an in/out parameter
inline int nextCode(int& x)
{
    if(x < 0)
    {
        return 0; // coding finished
    }
    else
    {
        int outCode;
        if(x & 1)
        {
            outCode = RUN_B;
        }
        else
        {
            outCode = RUN_A;
        }

        // Update X
        if(x < 2)
            x = -1; //signal encoding is finished
        else
            x = (x - 2) / 2;
        return outCode;
    }
}

};

#endif
