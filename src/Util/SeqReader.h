//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SeqReader - Reads fasta or fastq sequence files
//
#ifndef SEQREADER_H
#define SEQREADER_H

#include <fstream>
#include "Util.h"

enum RecordType
{
    RT_FASTA,
    RT_FASTQ,
    RT_UNKNOWN
};

enum SeqReaderFlag
{
    SRF_NO_VALIDATION,
    SRF_FULL_VALIDATION
};

//
class SeqReader
{
    public:
        SeqReader(std::string filename, SeqReaderFlag flag = SRF_FULL_VALIDATION);
        ~SeqReader();
        bool get(SeqRecord& sr);

    private:
        std::istream* m_pHandle;
        SeqReaderFlag m_flag;
};

#endif
