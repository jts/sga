//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLBWT - Run-length encoded Burrows Wheeler transform
//
#ifndef RLBWT_H
#define RLBWT_H

#include <deque>
#include "STCommon.h"
#include "Occurrence.h"
#include "SuffixArray.h"
#include "ReadTable.h"
#include "HitData.h"
#include "BWTReader.h"
#include "EncodedString.h"
#include "FMMarkers.h"
#include "RLUnit.h"
#include "HuffmanUtil.h"
#include "StreamEncoding.h"
#include "HuffmanForest.h"
#include "PackedTableDecoder.h"

// Defines
#define RLBWT_VALIDATE 1

typedef std::vector<uint8_t> RLRawData;

//
// RLBWT
//
class RLBWT
{
    public:
    
        // Constructors
        RLBWT(const std::string& filename, int sampleRate = DEFAULT_SAMPLE_RATE_SMALL);

        //    
        void setSampleRates(size_t largeSampleRate, size_t smallSampleRate);
        void initializeFMIndex(AlphaCount64& running_ac);

        // Append a symbol to the bw string
        void append(char b);

        inline char getChar(size_t idx) const
        {
            // Decompress stream up to the (idx + 1) character and return the last decompressed symbol
            size_t encoderIdx = 0;
            const LargeMarker marker = getLowerMarker(idx, encoderIdx);
            size_t current_position = marker.getActualPosition();
            size_t numToCount = idx - current_position + 1;
            //assert(numToCount < m_smallSampleRate);
            size_t symbol_index = marker.unitIndex;
            size_t numBitsRead = 0;

            char outBase = '\0';
            StreamEncode::SingleBaseDecode sbd(outBase);
            StreamEncode::decodeStream(&m_rlString[symbol_index], &m_rlString.back(), numToCount, numBitsRead, sbd);
            return outBase;
        }

        // Get the greatest interpolated marker whose position is less than or equal to position
        inline LargeMarker getLowerMarker(size_t position, size_t& encoderIdx) const
        {
            size_t target_small_idx = position >> m_smallShiftValue;
            if(target_small_idx >= m_smallMarkers.size())
                std::cout << "Pos: " << position << " targetsmall: " << target_small_idx << "\n";
            return getInterpolatedMarker(target_small_idx, encoderIdx);
        }

        // Return a LargeMarker with values that are interpolated by adding
        // the relative count nearest to the requested position to the last
        // LargeMarker
        inline LargeMarker getInterpolatedMarker(size_t target_small_idx, size_t& encoderIdx) const
        {
            // Calculate the position of the LargeMarker that the SmallMarker is relative to
            size_t target_position = target_small_idx << m_smallShiftValue;
            size_t curr_large_idx = target_position >> m_largeShiftValue;
            LargeMarker absoluteMarker = m_largeMarkers[curr_large_idx];
            assert(target_small_idx < m_smallMarkers.size());
            const SmallMarker& relative = m_smallMarkers[target_small_idx];
            alphacount_add16(absoluteMarker.counts, relative.counts);
            absoluteMarker.unitIndex += relative.unitCount;
            encoderIdx = relative.encoderIdx;
            return absoluteMarker;
        }

        inline BaseCount getPC(char b) const { return m_predCount.get(b); }

        // Return the number of times char b appears in bwt[0, idx]
        inline BaseCount getOcc(char b, size_t idx) const
        {
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            size_t encoderIdx = 0;
            const LargeMarker marker = getLowerMarker(idx, encoderIdx);
            size_t current_position = marker.getActualPosition();
            size_t numToCount = idx - current_position;
            assert(numToCount < m_smallSampleRate);
            size_t running_count = marker.counts.get(b);
            size_t symbol_index = marker.unitIndex;
            StreamEncode::BaseCountDecode bcd(b, running_count);
            size_t numBitsRead = 0;
            StreamEncode::decodeStream(&m_rlString[symbol_index], &m_rlString.back(), numToCount, numBitsRead, bcd);
            return running_count;
        }

        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline AlphaCount64 getFullOcc(size_t idx) const 
        { 
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            size_t encoderIdx = 0;
            const LargeMarker marker = getLowerMarker(idx, encoderIdx);
            size_t current_position = marker.getActualPosition();
            AlphaCount64 running_count = marker.counts;
            size_t numToCount = idx - current_position;

            assert(numToCount < m_smallSampleRate);
            size_t symbol_index = marker.unitIndex;
            StreamEncode::AlphaCountDecode acd(running_count);
            size_t numBitsRead = 0;
            StreamEncode::decodeStream(&m_rlString[symbol_index], &m_rlString.back(), numToCount, numBitsRead, acd);
            return running_count;
        }

        // Return the number of times each symbol in the alphabet appears ins bwt[idx0, idx1]
        inline AlphaCount64 getOccDiff(size_t idx0, size_t idx1) const 
        { 
            return getFullOcc(idx1) - getFullOcc(idx0); 
        }

        inline size_t getNumStrings() const { return m_numStrings; } 
        inline size_t getBWLen() const { return m_numSymbols; }
        inline size_t getNumRuns() const { return m_rlString.size(); }
        inline size_t getSmallSampleRate() const { return m_smallSampleRate; }
        // Return the first letter of the suffix starting at idx
        inline char getF(size_t idx) const
        {
            size_t ci = 0;
            while(ci < ALPHABET_SIZE && m_predCount.getByIdx(ci) <= idx)
                ci++;
            assert(ci != 0);
            return RANK_ALPHABET[ci - 1];
        }

        // Print the size of the BWT
        void printInfo() const;
        void print() const;
        void printRunLengths() const;
        
        void decodeToFile(const std::string& file);

        // IO
        friend class BWTReaderBinary;
        friend class BWTWriterBinary;
        friend class BWTReaderAscii;
        friend class BWTWriterAscii;

        // Default sample rates for the large (64-bit) and small (8-bit) occurrence markers
        static const int DEFAULT_SAMPLE_RATE_LARGE = 16384;
        static const int DEFAULT_SAMPLE_RATE_SMALL = 128;

    private:


        // Default constructor is not allowed
        RLBWT() {}
        
        // Calculate the number of markers to place
        size_t getNumRequiredMarkers(size_t n, size_t d) const;

        // The C(a) array
        AlphaCount64 m_predCount;
        
        // RL huffman tree
        HuffmanTreeCodec<int> m_rlHuffman;
        RLPackedTableDecoder m_rlDecodeTable;

        // The run-length encoded string
        RLRawData m_rlString;

        // The marker vectors
        LargeMarkerVector m_largeMarkers;
        SmallMarkerVector m_smallMarkers;

        // The number of strings in the collection
        size_t m_numStrings;

        // The total length of the bw string
        size_t m_numSymbols;

        // The sample rate used for the markers
        size_t m_largeSampleRate;
        size_t m_smallSampleRate;

        // The amount to shift values by to divide by m_sampleRate
        int m_smallShiftValue;
        int m_largeShiftValue;

};
#endif
