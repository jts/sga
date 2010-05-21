//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// EncodedString - Templated class to store
// a string from a reduced alphabet. The actual
// coding is done outside this class by the Codec
// passed in as the template parameter
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
            assert(false);
            // Deep copy
            _alloc(other.m_data, other.m_len);
        }

        //
        EncodedString(std::string seq)
        {
            _alloc(seq.c_str(), seq.length());
        }

        //
        ~EncodedString()
        {
            _dealloc();
        }

        // Operators
        EncodedString& operator=(const EncodedString& other)
        {
            assert(false);
            if(&other == this)
                return *this; // self-assign

            _dealloc();
            _alloc(other.m_data, other.m_len);
            return *this;
        }

        EncodedString& operator=(const std::string& str)
        {
            _dealloc();
            _alloc(str.c_str(), str.length());
            return *this;
        }

        size_t length() const
        {
            return m_len;
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

        // allocate storage and initialize the data
        void _alloc(const char* pData, size_t l)
        {
            m_len = l;
            // Get the number of units that need to be allocated from the codec
            size_t n_units = s_codec.getRequiredUnits(l);
            m_data = new StorageUnit[n_units](); // This zeros the memory
            m_capacity = s_codec.getCapacity(n_units);

            // Encode the string
            for(size_t i = 0; i < l; ++i)
                s_codec.store(m_data, i, pData[i]);
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
