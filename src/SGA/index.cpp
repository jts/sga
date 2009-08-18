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
	std::ifstream in(filename.c_str());
	ReadTable rt;

	std::string line;
	size_t count = 0;
	while(in >> line)
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
