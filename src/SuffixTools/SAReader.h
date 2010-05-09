//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SAReader.h - Read a suffix array file from disk
//
#ifndef SAREADER_H
#define SAREADER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"

const uint16_t SA_FILE_MAGIC = 0xCACA;

enum SAIOStage
{
	SAIOS_NONE,
	SAIOS_HEADER,
	SAIOS_ELEM,
	SAIOS_DONE
};

class SuffixArray;

class SAReader
{
	public:
		SAReader(const std::string& filename);
		~SAReader();

		//
		void read(SuffixArray* pSA);
		void readHeader(size_t& num_strings, size_t& num_elems);
		void readElems(SAElemVector& elemVector);
		SAElem readElem();

	private:
		std::istream* m_pReader;
		SAIOStage m_stage;
};

#endif
