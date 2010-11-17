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

// LargeMarker - To allow random access to the 
// BWT symbols and implement the occurrence array
// we keep a vector of symbol counts every D1 symbols.
// These counts are the absolute number of times each
// symbol has been seen up to that point.
// 
struct LargeMarker
{
    LargeMarker() : unitIndex(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getActualPosition() const
    {
        return counts.getSum();
    }

    void print() const
    {
        std::cout << "Large marker actual pos: " << getActualPosition() << "\n";
        std::cout << "Marker unit index: " << unitIndex << "\n";
        std::cout << "Marker counts: ";
        for(int i = 0; i < ALPHABET_SIZE; ++i)
        {
            std::cout << (int)counts.getByIdx(i) << " ";
        }
        std::cout << "\n";
    }    

    // Returns true if the data in the markers is identical
    bool operator==(const LargeMarker& rhs)
    {
        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(counts.getByIdx(i) != rhs.counts.getByIdx(i))
                return false;
        }
        return unitIndex == rhs.unitIndex;
    }

    // The number of times each symbol has been seen up to this marker
    AlphaCount64 counts; 

    // The index in the RLVector of the run that starts after
    // this marker. That is, if C = getActualPosition(), then
    // the run containing the B[C] is at unitIndex. This is not necessary
    // a valid index if there is a marker after the last symbol in the BWT
    size_t unitIndex;
};
typedef std::vector<LargeMarker> LargeMarkerVector;

// SmallMarker - Small markers contain the counts
// within an individual block of the BWT. In other words
// the small marker contains the count for the last D2 symbols
// 
struct SmallMarker
{
    SmallMarker() : unitCount(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getCountSum() const
    {
        return counts.getSum();
    }

    void print() const
    {
        for(int i = 0; i < ALPHABET_SIZE; ++i)
        {
            std::cout << (int)counts.getByIdx(i) << " ";
        }
        std::cout << "\n";
    }

    // The number of times each symbol has been seen up to this marker
    AlphaCount8 counts; 

    // The number of RL units in this block
    uint8_t unitCount;
};
typedef std::vector<SmallMarker> SmallMarkerVector;

//
// RLBWT
//
class RLBWT
{
    public:
    
        // Constructors
        RLBWT(const std::string& filename, int sampleRate = DEFAULT_SAMPLE_RATE);

        //    
        void initializeFMIndex();

        // Append a symbol to the bw string
        void append(char b);

        inline char getChar(size_t idx) const
        {
            // Calculate the Marker who's position is not less than idx
            const LargeMarker& upper = getUpperMarker(idx);
            size_t current_position = upper.getActualPosition();
            assert(current_position >= idx);

            size_t symbol_index = upper.unitIndex; 

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

        // Get the index of the marker nearest to position in the bwt
        inline size_t getNearestMarkerIdx(size_t position, size_t sampleRate, size_t shiftValue) const
        {
            size_t offset = MOD_POWER_2(position, sampleRate); // equivalent to position % sampleRate
            size_t baseIdx = position >> shiftValue;

            if(offset < (sampleRate >> 1))
            {
                return baseIdx;
            }
            else
            {
                return baseIdx + 1;
            }
        }        

        // Get the interpolated marker with position closest to position
        inline LargeMarker getNearestMarker(size_t position) const __attribute__((always_inline))

        {
            size_t nearest_small_idx = getNearestMarkerIdx(position, m_smallSampleRate, m_smallShiftValue);
            return getInterpolatedMarker(nearest_small_idx);
        }

        // Get the greatest interpolated marker whose position is less than or equal to position
        inline LargeMarker getLowerMarker(size_t position) const __attribute__((always_inline))
        {
            size_t target_small_idx = position >> m_smallShiftValue;
            return getInterpolatedMarker(target_small_idx);
        }

        // Get the lowest interpolated marker whose position is strictly greater than position
        inline LargeMarker getUpperMarker(size_t position) const __attribute__((always_inline))
        {
            size_t target_small_idx = (position >> m_smallShiftValue) + 1;
            return getInterpolatedMarker(target_small_idx);
        }

        // Return a LargeMarker with values that are interpolated by adding/subtracting all the SmallMarkers up to target_small_idx
        inline LargeMarker getInterpolatedMarker(size_t target_small_idx) const __attribute__((always_inline))
        {
            // Calculate the position of the LargeMarker closest to the target SmallMarker
            size_t target_position = target_small_idx << m_smallShiftValue;
            size_t curr_large_idx = getNearestMarkerIdx(target_position, m_largeSampleRate, m_largeShiftValue);
            size_t current_position = curr_large_idx << m_largeShiftValue;

            LargeMarker accumulated = m_largeMarkers[curr_large_idx];

            size_t current_small_idx = current_position >> m_smallShiftValue;

            // Each small block contains the count of symbols in the last sampleRate bases
            // When we count forward we skip the first block as it's count is included
            // in the LargeMarker at the same position. When counting backwards we include
            // the small marker at the same position of the LargeMarker which is why
            // the loop variable is incremented/decremented in different places in the loops below.
            if(current_position < target_position)
            {
                // Accumulate small blocks forward
                while(current_position < target_position)
                {
                    // Invariant: The accumulated count includes the small block ending
                    // at current_small_idx.
                    
                    // Add the counts in the next block 
                    current_small_idx += 1;
                    current_position += m_smallSampleRate;

                    const SmallMarker& small_marker = m_smallMarkers[current_small_idx];
                    // Add the small marker counts into the large marker
                    alphacount_add(accumulated.counts, small_marker.counts);
                    accumulated.unitIndex += small_marker.unitCount;
                }
            }
            else
            {
                // Accumulate small blocks backwards
                while(current_position > target_position)
                {
                    current_position -= m_smallSampleRate;
                    const SmallMarker& small_marker = m_smallMarkers[current_small_idx];
                    alphacount_subtract(accumulated.counts, small_marker.counts);
                    accumulated.unitIndex -= small_marker.unitCount;
                    current_small_idx -= 1;
                }
            }
            assert(current_position == target_position);
            return accumulated;
        }

        inline BaseCount getPC(char b) const { return m_predCount.get(b); }

        // Return the number of times char b appears in bwt[0, idx]
        inline BaseCount getOcc(char b, size_t idx) const
        {
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            const LargeMarker& marker = getNearestMarker(idx);
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
        inline AlphaCount64 getFullOcc(size_t idx) const 
        { 
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            const LargeMarker& marker = getNearestMarker(idx);
            size_t current_position = marker.getActualPosition();
            bool forwards = current_position < idx;

            AlphaCount64 running_count = marker.counts;
            size_t symbol_index = marker.unitIndex; 

            if(forwards)
                accumulateForwards(running_count, symbol_index, current_position, idx);
            else
                accumulateBackwards(running_count, symbol_index, current_position, idx);
            return running_count;
        }

        // Adds to the count of symbol b in the range [targetPosition, currentPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateBackwards(AlphaCount64& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
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
        inline void accumulateForwards(AlphaCount64& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const
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
        inline AlphaCount64 getOccDiff(size_t idx0, size_t idx1) const 
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
        void printRunLengths() const;

        // IO
        friend class BWTReaderBinary;
        friend class BWTWriterBinary;
        friend class BWTReaderAscii;
        friend class BWTWriterAscii;

        // Default sample rates for the large (64-bit) and small (8-bit) occurrence markers
        static const int DEFAULT_SAMPLE_RATE_LARGE = 1024;
        static const int DEFAULT_SAMPLE_RATE_SMALL = 128;
        static const int DEFAULT_SAMPLE_RATE = 128;

    private:


        // Default constructor is not allowed
        RLBWT() {}
        
        // Calculate the number of markers to place
        size_t getNumRequiredMarkers(size_t n, size_t d) const;

        // The C(a) array
        AlphaCount64 m_predCount;
        
        // The run-length encoded string
        RLVector m_rlString;

        // The marker vector
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
