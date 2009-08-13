#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include "Resolve.h"
#include "Util.h"

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string adjFile(argv[optind++]);
	
	ContigGraph* pContigGraph = createContigGraph(opt::k, contigFile, adjFile);

	//sg.stats();
	pContigGraph->validate();
	pContigGraph->writeDot("contigGraph.dot", DF_UNDIRECTED, colorVertexByAnnotation);
	//sg.simplify();
	//sg.stats();
	//sg.validate();
	delete pContigGraph;
}

// Default vertex color function, returns black for everything
std::string colorVertexByLength(Contig d)
{
	if(d.getLength() > 500)
	{
		return "red";
	}
	else
	{
		return "black";
	}
}

// Return the color based on the annotation of the contig
std::string colorVertexByAnnotation(Contig d)
{
	unsigned int annt = d.getAnnotation();
	
	int numSet = 0;
	std::string color = "black";

	if(annt & 0x1)
	{
		++numSet;
		color = "red";
	}
	if(annt & 0x2)
	{
		++numSet;
		color = "blue";
	}
	if(annt & 0x4)
	{
		++numSet;
		color = "green";
	}

	// Reference
	if(annt & 0x8)
	{
		++numSet;
		color = "black";
	}


	return color;
}





// 
// Handle command line arguments
//
void parseOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) 
	{
		std::cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 2) 
	{
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	} 
	else if (argc - optind > 3) 
	{
		std::cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) 
	{
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
}

