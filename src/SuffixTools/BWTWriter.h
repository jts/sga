//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter.h - Write a BWT file to disk
//
#ifndef BWTWRITER_H
#define BWTWRITER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTReader.h"

class BWT;

class BWTWriter
{
	public:
		BWTWriter(const std::string& filename);
		~BWTWriter();

		//
		void write(const BWT* pBWT);
		void writeHeader(const size_t& num_strings, const size_t& num_symbols);
		void writeBWStr(const std::string& str);
		void writeBWChar(char b);
		void writePred(const AlphaCount& pc);
		void writeOccurrence(const Occurrence& icc);

	private:
		std::ostream* m_pWriter;
		BWIOStage m_stage;
};

#endif
