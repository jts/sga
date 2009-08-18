//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// overlap - compute pairwise overlaps between reads
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"

int overlapMain(int argc, char** argv)
{
	if(argc == 0 || strcmp(argv[0], "help") == 0)
		printOverlapUsage();
	else
	{
		std::string indexFile(argv[0]);
		std::string readsFile(argv[1]);
		computeOverlaps(indexFile, readsFile);
		return 0;
	}
	return 0;
}

void computeOverlaps(std::string indexFile, std::string readsFile)
{
	std::ifstream inReads(readsFile.c_str());
	ReadTable rt;

	std::string line;
	size_t count = 0;
	while(inReads >> line)
	{
		if(count % 2 == 1)
		{
			Read r("", line);
			rt.addRead(r);
		}

		if(count % 10000 == 0)
			std::cout << "Processed " << count << "\n";
		++count;
	}
	std::cout << "Loaded " << rt.getCount() << " reads\n";

	// Load suffix array
	std::ifstream inSA(indexFile.c_str());
	SuffixArray sa;
	inSA >> sa;
	sa.validate(&rt);
}

void printOverlapUsage()
{
	std::cout << "sga overlap - compute pairwise overlaps between reads\n";
	std::cout << "Usage: sga overlap <index.sa> <in.fasta>\n";
}
