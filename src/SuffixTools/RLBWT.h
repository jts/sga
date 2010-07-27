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

#include "STCommon.h"
#include "Occurrence.h"
#include "SuffixArray.h"
#include "ReadTable.h"
#include "HitData.h"
#include "BWTReader.h"
#include "EncodedString.h"

// Defines
#define RL_COUNT_MASK 0x1F  //00011111
#define RL_SYMBOL_MASK 0xE0 //11100000
#define FULL_COUNT 31
#define RL_SYMBOL_SHIFT 5
#define RLBWT_VALIDATE 1

// A unit of the RLBWT is a pair of a symbol and its count
// The high 3 bits encodes the symbol to store
// The low 5 bits encodes the length of the run
struct RLUnit
{
    RLUnit() : data(0) {}
    RLUnit(char b) : data(1)
    {
        setChar(b);   
    }

    // Returns true if the count cannot be incremented
    inline bool isFull() const
    {
        return (data & RL_COUNT_MASK) == FULL_COUNT;
    }

    inline bool isEmpty() const
    {
        return (data & RL_COUNT_MASK) == 0;
    }

    inline bool isInitialized() const
    {
        return data > 0;
    }

    // 
    inline void incrementCount()
    {
#ifdef RLBWT_VALIDATE
        assert(!isFull());
#endif
        ++data;
    }

    // 
    inline void decrementCount()
    {
#ifdef RLBWT_VALIDATE
        assert(!isEmpty());
#endif
        --data;
    }    

    inline uint8_t getCount() const
    {
#ifdef RLBWT_VALIDATE
        assert((data & RL_COUNT_MASK) != 0);
#endif
        return data & RL_COUNT_MASK;
    }

    // Set the symbol
    inline void setChar(char symbol)
    {
        // Clear the current symbol
        data &= RL_COUNT_MASK;
        
        uint8_t code = BWT_ALPHABET::getRank(symbol);
        code <<= RL_SYMBOL_SHIFT;
        data |= code;
    }

    // Get the symbol
    inline char getChar() const
    {
        uint8_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return BWT_ALPHABET::getChar(code);
    }

    // 
    uint8_t data;

    friend class RLBWTReader;
    friend class RLBWTWriter;
};
typedef std::vector<RLUnit> RLVector;

// RLMarker - To allow random access to the 
// BWT symbols and implement the occurrence array
// we keep a vector of markers every D symbols.
// 
struct RLMarker
{
    RLMarker() : unitIndex(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getActualPosition() const
    {
        return counts.getSum();
    }

    // The number of times each symbol has been seen up to this marker
    AlphaCount counts; 

    // The index in the RLVector of the run that starts after
    // this marker. That is, if C = getActualPosition(), then
    // the run containing the B[C] is unitIndex. This is not necessary
    // a valid index if there is a marker after the last symbol in the BWT
    size_t unitIndex;
};
typedef std::vector<RLMarker> MarkerVector;

//
// RLBWT
//
class RLBWT
{
    public:
    
        // Constructors
        RLBWT(const std::string& filename, int sampleRate = DEFAULT_SAMPLE_RATE);
        RLBWT(const SuffixArray* pSA, const ReadTable* pRT);

        //    
        void initializeFMIndex();

        // Append a symbol to the bw string
        void append(char b);

        inline char getChar(size_t idx) const
        {
            // Calculate the Marker who's position is not less than idx
            const RLMarker& nearest = getUpperMarker(idx);
            size_t current_position = nearest.getActualPosition();
            assert(current_position >= idx);

            size_t symbol_index = nearest.unitIndex; 

            // Search backwards (towards 0) until idx is found
            while(current_position > idx)
            {
                assert(symbol_index != 0);
                symbol_index -= 1;
                current_position -= m_rlString[symbol_index].getCount();
            }

            // symbol_index is now the index of the run containing the idx symbol
            const RLUnit& unit = m_rlString[symbol_index];
            assert(current_position <= idx && current_position + unit.getCount() >= idx);
            return unit.getChar();
        }
    
        // Get the first marker with a position that is guaranteed to be 
        // no greater than idx
        inline const RLMarker& getUpperMarker(size_t idx) const
        {
            //printf("UPPER idx: %zu shifted: %zu nm: %zu\n", idx, idx >> m_shiftValue, m_markers.size());
            idx >>= m_shiftValue;
            ++idx;
#ifdef RLBWT_VALIDATE
            assert(idx < m_markers.size());
#endif
            return m_markers[idx];
        }

        // Get the the marker who's position is estimated to be at idx
        // but it may be slightly more
        inline const RLMarker& getLowerMarker(size_t idx) const
        {
            //printf("LOWER idx: %zu shifted: %zu nm: %zu\n", idx, idx >> m_shiftValue, m_markers.size());
            idx >>= m_shiftValue;
#ifdef RLBWT_VALIDATE
            assert(idx < m_markers.size());
#endif
            return m_markers[idx];
        }

        // Get the nearest marker to idx.
        inline const RLMarker& getNearestMarker(size_t idx) const
        {
            size_t offset = MOD_POWER_2(idx, m_sampleRate); // equivalent to idx % m_sampleRate
            if(offset < (m_sampleRate >> 1))
            {
                // Choose lower marker
                return getLowerMarker(idx);    
            }
            else
            {
                // Choose upper marker
                return getUpperMarker(idx);
            }
        }

        inline BaseCount getPC(char b) const { return m_predCount.get(b); }

        // Return the number of times char b appears in bwt[0, idx]
        inline BaseCount getOcc(char b, size_t idx) const
        {
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            const RLMarker& marker = getNearestMarker(idx);
            size_t current_position = marker.getActualPosition();
            bool forwards = current_position < idx;
            //printf("cp: %zu idx: %zu f: %d dist: %d\n", current_position, idx, forwards, (int)idx - (int)current_position);

            size_t running_count = marker.counts.get(b);
            size_t symbol_index = marker.unitIndex; 

            if(forwards)
                accumulateForwards(b, running_count, symbol_index, current_position, idx);
            else
                accumulateBackwards(b, running_count, symbol_index, current_position, idx);
            return running_count;
        }

        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline AlphaCount getFullOcc(size_t idx) const 
        { 
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            const RLMarker& marker = getNearestMarker(idx);
            size_t current_position = marker.getActualPosition();
            bool forwards = current_position < idx;

            AlphaCount running_count = marker.counts;
            size_t symbol_index = marker.unitIndex; 

            if(forwards)
                accumulateForwards(running_count, symbol_index, current_position, idx);
            else
                accumulateBackwards(running_count, symbol_index, current_position, idx);
            return running_count;
        }

        // Adds to the count of symbol b in the range [targetPosition, currentPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateBackwards(AlphaCount& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
        {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition)
            {
                size_t diff = currentPosition - targetPosition;
#ifdef RLBWT_VALIDATE                
                assert(currentUnitIndex != 0);
#endif
                --currentUnitIndex;

                const RLUnit& curr_unit = m_rlString[currentUnitIndex];
                uint8_t curr_count = curr_unit.getCount();
                if(curr_count > diff)
                    curr_count = diff;
                
                char curr_base = curr_unit.getChar();
                running_count.subtract(curr_base, curr_count);
                currentPosition -= curr_count;
            }
        }

        // Adds to the count of symbol b in the range [currentPosition, targetPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateForwards(AlphaCount& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
        {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition)
            {
                size_t diff = targetPosition - currentPosition;
#ifdef RLBWT_VALIDATE
                assert(currentUnitIndex != m_rlString.size());
#endif
                const RLUnit& curr_unit = m_rlString[currentUnitIndex];
                uint8_t curr_count = curr_unit.getCount();
                if(curr_count > diff)
                    curr_count = diff;
                
                char curr_base = curr_unit.getChar();
                running_count.add(curr_base, curr_count);
                currentPosition += curr_count;
                ++currentUnitIndex;
            }
        }

        // Adds to the count of symbol b in the range [targetPosition, currentPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateBackwards(char b, size_t& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
        {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition)
            {
                size_t diff = currentPosition - targetPosition;
#ifdef RLBWT_VALIDATE                
                assert(currentUnitIndex != 0);
#endif
                --currentUnitIndex;

                const RLUnit& curr_unit = m_rlString[currentUnitIndex];
                uint8_t curr_count = curr_unit.getCount();
                if(curr_count > diff)
                    curr_count = diff;
                
                char curr_base = curr_unit.getChar();
                if(curr_base == b)
                    running_count -= curr_count;
                currentPosition -= curr_count;
            }
        }

        // Adds to the count of symbol b in the range [currentPosition, targetPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateForwards(char b, size_t& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
        {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition)
            {
                size_t diff = targetPosition - currentPosition;
#ifdef RLBWT_VALIDATE
                assert(currentUnitIndex != m_rlString.size());
#endif
                const RLUnit& curr_unit = m_rlString[currentUnitIndex];
                uint8_t curr_count = curr_unit.getCount();
                if(curr_count > diff)
                    curr_count = diff;
                
                char curr_base = curr_unit.getChar();
                if(curr_base == b)
                    running_count += curr_count;
                currentPosition += curr_count;
                ++currentUnitIndex;
            }
        }

        // Return the number of times each symbol in the alphabet appears ins bwt[idx0, idx1]
        inline AlphaCount getOccDiff(size_t idx0, size_t idx1) const 
        { 
            return getFullOcc(idx1) - getFullOcc(idx0); 
        }

        inline size_t getNumStrings() const { return m_numStrings; } 
        inline size_t getBWLen() const { return m_numSymbols; }
        inline size_t getNumRuns() const { return m_rlString.size(); }

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
        void validate() const;

        // IO
        friend class BWTReaderBinary;
        friend class BWTWriterBinary;
        friend class BWTReaderAscii;
        friend class BWTWriterAscii;

        static const int DEFAULT_SAMPLE_RATE = 128;

    private:


        // Default constructor is not allowed
        RLBWT() {}

        // The C(a) array
        AlphaCount m_predCount;
        
        // The run-length encoded string
        RLVector m_rlString;

        // The marker vector
        MarkerVector m_markers;

        // The number of strings in the collection
        size_t m_numStrings;

        // The total length of the bw string
        size_t m_numSymbols;

        // The sample rate used for the markers
        size_t m_sampleRate;

        // The amount to shift values by to divide by m_sampleRate
        int m_shiftValue;

};
#endif
