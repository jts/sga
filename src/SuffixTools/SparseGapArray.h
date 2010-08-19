//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SparseGapArray
//
// The sparse gap array uses a small integer type to
// keep the count for each element. A small number
// of positions will have a count that is larger
// than the maximum value the integer can represent.
// These values are automatically overflowed into a 
// hash map with a larger integer type. In the typical
// case this will require far less memory than the
// SimpleGapArray. The storage unit is templated
// to allow a dynamic tradeoff between the amount
// of memory required to store the simple cases 
// and the overflow rate.
//
#ifndef SPARSEGAPARRAY_H
#define SPARSEGAPARRAY_H

#include "HashMap.h"
#include "GapArray.h"

// Template base storage for the sparse gap array
template<class IntType>
class SparseBaseStorage
{
    public:
        //
        void resize(size_t n)
        {
            m_data.resize(n);
        }

        //
        inline void set(size_t i, IntType c)
        {
            m_data[i] = c;
        }

        //
        inline void increment(size_t i)
        {
            ++m_data[i];
        }

        //
        inline IntType get(size_t i) const
        {
            return m_data[i];
        }

        inline static size_t getMax()
        {
            static size_t max = std::numeric_limits<IntType>::max();
            return max;
        }

        //
        size_t size() const
        {
            return m_data.size();
        }

    private:
        //
        typedef std::vector<IntType> BaseStorageVector;
        BaseStorageVector m_data;
};

static const uint8_t storage4_offset_mask[]={0xF0,0x0F}; // 11110000, 00001111
static const uint8_t storage4_shift[]={4, 0};

// specialized 4-bit per value base storage for the sparse gap array
class SparseBaseStorage4
{
    public:

        SparseBaseStorage4() : m_numElems(0) {}

        //
        void resize(size_t n)
        {
            // With 4 bits per entry, we need n / 2 bytes of storage
            size_t numBytes = ((n & 1) == 0) ? n / 2 : (n / 2) + 1;
            m_data.resize(numBytes);
            m_numElems = n;
        }

        //
        inline void set(size_t i, uint8_t c)
        {
            assert(i < m_numElems);
            assert(c <= getMax());

            size_t idx = i >> 1; // equivalent to i / 2
            size_t offset = i & 1; // 0 if even, 1 if odd
            uint8_t mask = storage4_offset_mask[offset];
            uint8_t& elem = m_data[idx];

            // Clear current value
            elem &= ~mask;
            
            // Set new value
            elem |= c << storage4_shift[offset];
        }

        //
        inline uint8_t get(size_t i) const
        {
            assert(i < m_numElems);

            size_t idx = i >> 1; // equivalent to i / 2
            size_t offset = i & 1; // 0 if even, 1 if odd
            uint8_t mask = storage4_offset_mask[offset];
            uint8_t elem = m_data[idx];

            elem &= mask;
            elem >>= storage4_shift[offset];
            assert(elem <= getMax());
            return elem;
        }

        inline static size_t getMax()
        {
            return 15;
        }

        //
        size_t size() const
        {
            return m_numElems;
        }

    private:
        
        //

        typedef std::vector<uint8_t> BaseStorageVector;
        BaseStorageVector m_data;
        size_t m_numElems;
};

template<class BaseStorage, class OverflowStorage>
class SparseGapArray : public GapArray
{
    public:
        SparseGapArray() {}
        ~SparseGapArray()
        {
            //printf("SparseGapArray -- n: %zu overflow: %zu (%lf)\n", m_baseStorage.size(), m_overflow.size(), (double)m_overflow.size() / m_baseStorage.size());    
        }

        //
        void resize(size_t n)
        {
            m_baseStorage.resize(n);
        }

        //
        void increment(size_t i)
        {
            assert(i < m_baseStorage.size());    
            size_t count = m_baseStorage.get(i);
            if(count == getBaseMax())
            {
                // Increment overflow
                incrementOverflow(i);
            }
            else
            {
                // Increment the base count
                // If it is now the maximum representable value
                // add it to the overlow map
                ++count;
                m_baseStorage.set(i, count);

                if(count == getBaseMax())
                {
                    initOverflow(i, count);
                }
            }
        }

        //
        void initOverflow(size_t i, OverflowStorage c)
        {
            m_overflow.insert(std::make_pair<size_t, OverflowStorage>(i, c));
        }

        //
        void incrementOverflow(size_t i)
        {
            // precondition: there must be an entry in the hash for this key
            assert(m_overflow.count(i) != 0);
            ++m_overflow[i];
        }

        //
        size_t get(size_t i) const
        {
            size_t count = m_baseStorage.get(i);
            if(count == getBaseMax())
            {
                typename OverflowHash::const_iterator iter = m_overflow.find(i);
                assert(iter != m_overflow.end());
                return iter->second;
            }
            else
            {
                return count;
            }
        }

        //
        size_t getBaseMax() const
        {
            return BaseStorage::getMax();
        }

        //
        size_t size() const
        {
            return m_baseStorage.size();
        }

   private:

        typedef SparseHashMap<size_t, OverflowStorage> OverflowHash;
        OverflowHash m_overflow;
        BaseStorage m_baseStorage;
};

typedef SparseGapArray<SparseBaseStorage4, size_t> SparseGapArray4;
typedef SparseGapArray<SparseBaseStorage<uint8_t>, size_t> SparseGapArray8;
typedef SparseGapArray<SparseBaseStorage<uint16_t>, size_t> SparseGapArray16;
typedef SparseGapArray<SparseBaseStorage<uint32_t>, size_t> SparseGapArray32;

#endif
