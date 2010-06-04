//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter.h - Write a BWT file to disk
//
#ifndef BWTWRITER_H
#define BWTWRITER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTReader.h"
#include "EncodedString.h"

class SBWT;
class RLBWT;

class BWTWriter
{
    public:
        BWTWriter(const std::string& filename);
        ~BWTWriter();

        //
        void write(const SBWT* pBWT);
        void write(const RLBWT* pRLBWT);

        void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag);
        void writeBWStr(const BWTString& str); 
        void writeBWChar(char b);
        void writePred(const AlphaCount& pc);
        void writeOccurrence(const Occurrence& icc);

    private:
        std::ostream* m_pWriter;
        BWIOStage m_stage;
};

#endif
