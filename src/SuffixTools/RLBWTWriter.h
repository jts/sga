//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLBWTWriter - Write a run-length encoded BWT to disk
//
#ifndef RLBWTWRITER_H
#define RLBWTWRITER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTReader.h"
#include "RLBWTReader.h"
#include "EncodedString.h"
#include "RLBWT.h"

class SBWT;
class RLBWT;

class RLBWTWriter
{
    public:
        RLBWTWriter(const std::string& filename);
        ~RLBWTWriter();

        // Write an RLBWT file directly from a suffix array and read table
        void write(const SuffixArray* pSA, const ReadTable* pRT);
        void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag);
        void writeBWChar(char b);
        void finalize(); // this method must be called after writing the BW string

    private:

        void writeRun(RLUnit& unit);

        std::ostream* m_pWriter;
        size_t m_numRuns;
        std::streampos m_runFileOffset;
        RLUnit m_currRun;
        BWIOStage m_stage;
};

#endif
