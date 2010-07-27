//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLRLBWTReader - Read a run length encoded BWT file from disk
//
#ifndef RLBWTREADER_H
#define RLBWTREADER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "EncodedString.h"
#include "BWTReader.h"
#include "RLBWT.h"

const uint16_t RLBWT_FILE_MAGIC = 0xCACA;

class SBWT;
class RLBWT;

class RLBWTReader
{
    public:
        RLBWTReader(const std::string& filename);
        ~RLBWTReader();

        //
        void read(RLBWT* pRLBWT);
        void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        char readBWChar();
        void readRuns(RLVector& out, size_t numRuns);

    private:
        std::istream* m_pReader;
        BWIOStage m_stage;
        RLUnit m_currRun;
        size_t m_numRunsOnDisk;
        size_t m_numRunsRead;
};

#endif
