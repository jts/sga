//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Quality - functions for manipulating quality values
//
#ifndef QUALITY_H
#define QUALITY_H

#include <inttypes.h>
#include <vector>
#include <math.h>
#include <assert.h>
#include <string>
#include <iostream>

static const int DEFAULT_QUAL_SCORE = 15;

typedef std::vector<double> DoubleVector;
namespace Quality
{
    // Phred score transformations
    inline int char2phred(char b)
    {
        uint8_t v = b;
        assert(v >= 33);
        return v - 33;
    }

    inline char phred2char(int p)
    {
        uint8_t v = (p <= 93) ? p : 93;
        char c = v + 33;
        return c;
    }

    inline int prob2phred(double p)
    {
        return static_cast<int>(round(-10.0f * log10(p)));
    }

    inline int lnprob2phred(double lp)
    {
        // change base
        static double transform = log(10);
        lp /= transform;
        return static_cast<int>(round(-10.0f * lp));
    }

    // Return a uniform log-scaled quality vector of the given size
    DoubleVector uniformLogProbVector(double p_error, size_t n);

    // Return a string which encodes the log-scaled double vector as characters
    std::string encodeLogProbVector(const DoubleVector& dv);
};


#endif
