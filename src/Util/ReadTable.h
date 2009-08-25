//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#ifndef READTABLE_H
#define READTABLE_H
#include "Util.h"
#include "SeqReader.h"

typedef std::vector<SeqItem> ReadVector;

class ReadTable
{
	public:
		//
		ReadTable() {}
		ReadTable(std::string filename);
		
		//
		void initializeReverse(const ReadTable* pRT);
		void addRead(const SeqItem& r);
		
		//
		const SeqItem& getRead(size_t idx) const;
		size_t getReadLength(size_t idx) const;
		size_t getCount() const;
		size_t getSumLengths() const;

		//
		friend std::ostream& operator<<(std::ostream& out, const ReadTable& rt);


	private:
		ReadVector m_table;
};

#endif
