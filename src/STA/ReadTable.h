//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#ifndef READDICTIONARY_H
#define READDICTIONARY_H
#include "STCommon.h"

typedef std::string DNASequence;

struct Read
{
	Read(std::string s1, std::string s2) : id(s1), seq(s2) {}
	std::string id;
	DNASequence seq;
};

typedef std::vector<Read> ReadVector;

class ReadTable
{
	public:
		ReadTable() {}
		
		void addRead(const Read& r);
		const Read& getRead(size_t idx) const;
		size_t getCount() const;

	private:
		ReadVector m_table;
};

#endif
