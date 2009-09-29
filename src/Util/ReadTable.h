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
#include <map>

typedef std::vector<SeqItem> ReadVector;
typedef std::map<std::string, SeqItem*> ReadIndex;

class ReadTable
{
	public:
		//
		ReadTable() : m_pIndex(NULL) {}
		ReadTable(std::string filename);
		~ReadTable();

		//
		void initializeReverse(const ReadTable* pRT);
		void addRead(const SeqItem& r);
		void indexReadsByID();

		//
		const SeqItem& getRead(size_t idx) const;
		const SeqItem& getRead(const std::string& id) const;
		size_t getReadLength(size_t idx) const;
		size_t getCount() const;
		size_t getSumLengths() const;

		//
		friend std::ostream& operator<<(std::ostream& out, const ReadTable& rt);


	private:
		ReadVector m_table;
		// Index of readid -> SeqItem
		// It is not build be default to save memory
		ReadIndex* m_pIndex; 
};

#endif
