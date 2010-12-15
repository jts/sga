//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTCompressor - Compress a BWT to disk
//
#ifndef BWTCOMPRESSOR_H
#define BWTCOMPRESSOR_H

#include <deque>
#include <iostream>
#include <map>
#include "HuffmanTreeCodec.h"
#include "Alphabet.h"
#include "FMMarkers.h"

// typedefs
typedef std::deque<char> CharDeque;
typedef std::map<char, int> CharIntMap;

//
class BWTCompressor
{
    public:
        
        BWTCompressor();
        ~BWTCompressor();

        void initialize(const std::string& filename, size_t largeSampleRate, 
                        size_t smallSampleRate, size_t numSymbols);

        void writeSymbol(char symbol, std::ostream* pWriter);
        void flush(std::ostream* pWriter);
        void writeLargeMarkers(std::ostream* pWriter);
        void writeSmallMarkers(std::ostream* pWriter);

        size_t getNumBytesWrote() const;
        size_t getNumSymbolsWrote() const;

        const LargeMarkerVector& getLargeMarkerVector() const;
        const SmallMarkerVector& getSmallMarkerVector() const;
        AlphaCount64 getRunningCount() const;
        
    private:

        void buildMarkers();
        void binaryFileCopy(std::istream* pReader, std::ostream* pWriter, 
                            size_t numBytes, size_t unitSize);

        size_t m_largeSampleRate;
        size_t m_smallSampleRate;
        size_t m_expectedSymbols;

        // 
        size_t m_symbolsWrote;
        size_t m_unitsWrote;
        AlphaCount64 m_runningAC;
        CharDeque m_symbolBuffer;

        //
        HuffmanTreeCodec<int> m_rlEncoder;

        // The Large/Small markers are written to temporary
        // files then copied to the bwt file
        std::string m_tempLargeFilename;
        std::string m_tempSmallFilename;

        std::ostream* m_pTempLargeMarkerWriter;
        std::ostream* m_pTempSmallMarkerWriter;

        // The last large/small marker that were written out
        SmallMarker m_prevSmallMarker;
        LargeMarker m_prevLargeMarker;

        size_t m_numSmallMarkersWrote;
        size_t m_numLargeMarkersWrote;
};

#endif
