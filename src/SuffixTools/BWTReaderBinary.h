//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Read a run length encoded binary BWT file from disk
//
#ifndef BWTREADERBINARY_H
#define BWTREADERBINARY_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "EncodedString.h"
#include "BWTReader.h"
#include "RLBWT.h"

class SBWT;
class RLBWT;

class BWTReaderBinary : public IBWTReader
{
    public:
        BWTReaderBinary(const std::string& filename);
        virtual ~BWTReaderBinary();

        //
        virtual void read(RLBWT* pRLBWT);
        virtual void read(SBWT* pSBWT);

        virtual void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        virtual char readBWChar();

    private:
    
        void readCompressedData(RLBWT* pBWT, size_t numUnits);
        void readMarkers(RLBWT* pBWT);
        void readPredCounts(RLBWT* pBWT);

        // 
        void decompressBlock();
        SmallMarker readSmallMarker();
        void fillCompressedBuffer();
        void _decompressBuffer(const SmallMarker& smallMarker);

        std::istream* m_pReader;
        BWIOStage m_stage;
        RLUnit m_currRun;
        size_t m_numUnitsOnDisk;
        size_t m_smallSampleRate;
        size_t m_largeSampleRate;

        // Data members required for streaming one symbol 
        // at a time.
        size_t m_numUnitsRead;
        size_t m_numSymbolsTotal;

        size_t m_totalSmallMarkers;
        size_t m_readSmallMarkers;
        size_t m_smallMarkerOffset;
        size_t m_decompressedTotal;

        RLPackedTableDecoder m_rlDecodeTable;
        CharDeque m_decompressedBuffer;
        std::vector<uint8_t> m_compressedBuffer;
};

#endif
