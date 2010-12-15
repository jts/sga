//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Write a run-length encoded BWT to a binary file
//
#ifndef BWTWRITERBINARY_H
#define BWTWRITERBINARY_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTWriter.h"
#include "BWTReaderBinary.h"
#include "EncodedString.h"
#include "RLBWT.h"
#include "BWTCompressor.h"

class SBWT;
class RLBWT;

class BWTWriterBinary : public IBWTWriter
{
    public:
        BWTWriterBinary(const std::string& filename, int smallSampleRate);
        virtual ~BWTWriterBinary();

        // Write an RLBWT file directly from a suffix array and read table
        virtual void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag);
        virtual void writeBWChar(char b);
        virtual void finalize(); // this method must be called after writing the BW string

    private:

        std::ostream* m_pWriter;
        size_t m_smallSampleRate;
        std::string m_filename;
        size_t m_numRuns;
        std::streampos m_headerFileOffset;
        std::streampos m_largeMarkerOffset;
        std::streampos m_smallMarkerOffset;
        RLUnit m_currRun;
        BWIOStage m_stage;
        BWTCompressor m_compressor;
};

#endif
