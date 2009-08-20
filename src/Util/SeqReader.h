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

//
class SeqReader
{
	public:
		SeqReader(std::string filename);
		~SeqReader();
		bool get(SeqItem& si);

	private:
		std::ifstream m_fileHandle;
};

#endif
