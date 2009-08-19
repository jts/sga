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
#include "BWT.h"

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
	ReadTable rt(readsFile);

	// Load suffix array
	std::ifstream inSA(indexFile.c_str());
	SuffixArray sa;
	inSA >> sa;
	sa.validate(&rt);

	// Convert SA to a BWT
	BWT b(&sa, &rt);
	b.print(&rt);

	// Compute overlaps
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		b.getOverlaps(rt.getRead(i).seq, 30);
	}
}

void printOverlapUsage()
{
	std::cout << "sga overlap - compute pairwise overlaps between reads\n";
	std::cout << "Usage: sga overlap <index.sa> <in.fasta>\n";
}
