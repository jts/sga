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

        void initialize(size_t largeSampleRate, size_t smallSampleRate, size_t numSymbols);

        void writeSymbol(char symbol, std::ostream* pWriter);
        void flush(std::ostream* pWriter);
        void writeMarkers();

        size_t getNumBytesWrote() const;
        const LargeMarkerVector& getLargeMarkerVector() const;
        const SmallMarkerVector& getSmallMarkerVector() const;
        
    private:

        // Append markers to the stream and return
        // the last small marker
        SmallMarker& appendMarkers();

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

        LargeMarkerVector m_largeMarkers;
        SmallMarkerVector m_smallMarkers;
};

#endif
