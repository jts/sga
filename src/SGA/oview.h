//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// oview - view overlaps between reads
//
#ifndef OVIEW_H
#define OVIEW_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"

// typedefs
typedef std::map<std::string, OverlapVector> OverlapMap;

struct DrawData
{
	DrawData(int o, std::string n, std::string s, int nd, int ol) : offset(o), name(n), seq(s), numDiff(nd), overlapLen(ol) {}
	int offset;
	std::string name;
	std::string seq;
	int numDiff;
	int overlapLen;

	bool operator<(const DrawData& other) const
	{
		return offset < other.offset;
	}
};

typedef std::vector<DrawData> DrawVector;

// functions

int oviewMain(int argc, char** argv);
void correctReads(const ReadTable* pRT, const OverlapMap* pOM);
std::string correct(const SeqItem& read, const ReadTable* pRT, const OverlapMap* pOM, int& num_corrected);
void drawAlignment(std::string rootID, const ReadTable* pRT, const OverlapMap* pOM);
void drawMulti(std::string rootName, int root_len, DrawVector& dv);
void parseOverlaps(std::string filename, OverlapMap& overlapMap);
void parseOviewOptions(int argc, char** argv);

#endif
