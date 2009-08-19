//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// index - index reads using a suffix tree
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "index.h"
#include "SuffixArray.h"
#include "SeqReader.h"

int indexMain(int argc, char** argv)
{
	if(argc != 1)
	{
		printIndexUsage();
		return 1;
	}
	else
	{
		std::string filename(argv[0]);
		std::cout << "Building index for " << filename << "\n";
		buildIndex(filename);
		return 0;
	}
}

void buildIndex(std::string filename)
{
	ReadTable rt(filename);

	// Make initial suffix arrays
	SuffixArray sa;
	sa.initialize(rt);
	sa.sort(&rt);
	sa.validate(&rt);

	std::ofstream out("indexed.sa");
	out << sa;
	out.close();
}

void printIndexUsage()
{
	std::cout << "sga index - create a suffix array for the given reads file\n";
	std::cout << "Usage: sga index <in.fasta>\n";
}
