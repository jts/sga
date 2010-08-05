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

/* 
// This is a simple wrapper around one of the standard 
// integer representations, uint8_t, uint16_t, etc
template<class IntType>
class SimpleIntStorage
{
    IntType() : m_value(0) {}
    
    void increment() { ++m_value; }
    IntType get() { return m_value; }
    void set(IntType v) { m_value = v; }
    
    // Return the maximum representable value with IntType
    static size_t getMax() const
    {
        static size_t max = std::numeric_limits<IntType>::max();
        return max;
    }
    
    bool isMax() const 
    { 
        if(m_value >= getMax())
            return true;
        else
            return false;
    }

    IntType m_value;
};
*/

template<class BaseStorage, class OverflowStorage>
class SparseGapArray : public GapArray
{
    public:
        SparseGapArray() {}
        ~SparseGapArray()
        {
            printf("SparseGapArray -- n: %zu overflow: %zu (%lf)\n", m_baseCounts.size(), m_overflow.size(), (double)m_overflow.size() / m_baseCounts.size());    
        }

        //
        void resize(size_t n)
        {
            m_baseCounts.resize(n);
        }

        //
        void increment(size_t i)
        {
            assert(i < m_baseCounts.size());    
            size_t count = m_baseCounts[i];
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
                m_baseCounts[i] = count;

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
            size_t count = m_baseCounts[i];
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
            return std::numeric_limits<BaseStorage>::max();
        }

        //
        size_t size() const
        {
            return m_baseCounts.size();
        }

   private:
        typedef std::vector<BaseStorage> BaseVector;
        typedef SparseHashMap<size_t, OverflowStorage> OverflowHash;

        OverflowHash m_overflow;
        BaseVector m_baseCounts;
};

typedef SparseGapArray<uint8_t, size_t> SparseGapArray8;

#endif
