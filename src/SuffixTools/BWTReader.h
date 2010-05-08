//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReader.h - Read a BWT file from disk
//
#ifndef BWTREADER_H
#define BWTREADER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"

const uint16_t BWT_FILE_MAGIC = 0xEFEF;

enum BWIOStage
{
	IOS_NONE,
	IOS_HEADER,
	IOS_BWSTR,
	IOS_PC,
	IOS_OCC,
	IOS_DONE
};

enum BWFlag
{
	BWF_NOFMI = 0,
	BWF_HASFMI
};

class BWT;

class BWTReader
{
	public:
		BWTReader(const std::string& filename);
		~BWTReader();

		//
		void read(BWT* pBWT);
		void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
		void readBWStr(std::string& out_str);
		char readBWChar();
		void readPred(AlphaCount& out_pc);
		void readOccurrence(Occurrence& out_icc);

	private:
		std::istream* m_pReader;
		BWIOStage m_stage;
};

#endif
