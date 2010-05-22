//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// EncodedString - Templated class to store
// a string from a reduced alphabet. The actual
// encoding is done outside this class by the Codec
// template parameter
//
#ifndef ENCODEDSTRING_H
#define ENCODEDSTRING_H
#include <string.h>
#include <string>
#include <assert.h> 
#include "DNACodec.h"

template<class Codec>
class EncodedString
{
        typedef typename Codec::UNIT_TYPE StorageUnit;

    public:
        
        // Constructors/Destructors
        EncodedString() : m_len(0), m_data(0) {}

        //
        EncodedString(const EncodedString& other)
        {
            // deep copy
            _alloc(other.m_len);
            _copy(other);
        }

        //
        EncodedString(const std::string& seq)
        {
            size_t n = seq.length();
            _alloc(n);
            _copy(seq.c_str(), n);
        }

        //
        ~EncodedString()
        {
            _dealloc();
        }

        // Operators
        EncodedString& operator=(const EncodedString& other)
        {
            if(&other == this)
                return *this; // self-assign

            _dealloc();
            _alloc(other.m_len);
            _copy(other);
            return *this;
        }

        EncodedString& operator=(const std::string& str)
        {
            size_t n = str.length();
            _dealloc();
            _alloc(n);
            _copy(str.c_str(), n);
            return *this;
        }

        size_t length() const
        {
            return m_len;
        }

        size_t capacity() const
        {
            return m_capacity;
        }

        bool empty() const
        {
            return m_len == 0;
        }

        // Get the character at the given position
        inline char get(size_t idx) const
        {
            return s_codec.get(m_data, idx);
        }

        //
        std::string toString() const
        {
            std::string out(m_len, 'A');
            for(size_t i = 0; i < m_len; ++i)
                out[i] = s_codec.get(m_data, i);
            return out;
        }

        //
        friend std::ostream& operator<<(std::ostream& out, const EncodedString<Codec>& a)
        {
            out << "Length: " << a.m_len << "\n";
            out << "Capacity: " << a.m_capacity << "\n";
            out << "Str: " << a.toString() << "\n";
            return out;
        }

    private:

        // functions
        void _copy(const char* pData, size_t n)
        {
            // this assumes that storage for n characters has been 
            // allocated
            assert(m_capacity >= n);

            for(size_t i = 0; i < n; ++i)
                s_codec.store(m_data, i, pData[i]);
        }
        
        //
        void _copy(const EncodedString& other)
        {
            // storage should have been allocated already
            assert(m_capacity = other.m_capacity);
            size_t n = s_codec.getRequiredUnits(other.m_capacity);
            for(size_t i = 0; i < n; ++i)
                m_data[i] = other.m_data[i];
        }

        // allocate storage for n symbols
        void _alloc(size_t n)
        {
            m_len = n;
            // Get the number of units that need to be allocated from the codec
            size_t n_units = s_codec.getRequiredUnits(n);
            m_data = new StorageUnit[n_units](); // This zeros the memory
            m_capacity = s_codec.getCapacity(n_units);
        }

        // deallocate storage
        void _dealloc()
        {
            if(m_data != NULL)
                delete [] m_data;
            m_data = 0;
            m_len = 0;
            m_capacity = 0;
        }

        // data
        static Codec s_codec;
        size_t m_len;
        size_t m_capacity;
        StorageUnit* m_data;
};

#endif
