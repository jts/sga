//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadInfoTable - A 0-indexed table of ID, length pairs
// Used to convert suffix array hits to overlaps
//
#ifndef READINFOTABLE_H
#define READINFOTABLE_H
#include "Util.h"
#include "SeqReader.h"
#include <map>

struct ReadInfo
{
    ReadInfo(const std::string& i, const uint32_t& l) : id(i), length(l) {}
    std::string id;
    uint32_t length;
};

typedef std::vector<ReadInfo> InfoVector;

class ReadInfoTable
{
    public:
        //
        ReadInfoTable() {}

        // Load the table using the read in filename
        // If num_expected > 0, reserve room in the table for num_expected reads
        ReadInfoTable(std::string filename, size_t num_expected = 0);
        ~ReadInfoTable();

        //
        void addReadInfo(const ReadInfo& r);
        const ReadInfo& getReadInfo(size_t idx) const;
        std::string getReadID(size_t idx) const;
        size_t getReadLength(size_t idx) const;
        size_t getCount() const;
        size_t countSumLengths() const;
        void clear();

    private:
        InfoVector m_table;
};

#endif
