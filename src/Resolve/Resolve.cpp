#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include "Resolve.h"
#include "SeqGraph.h"
#include "SeqVertex.h"
#include "Util.h"
#include "Contig.h"

struct AdjInfo
{
	AdjInfo(int o) : overlap(o) {}
	int overlap;
	Edge e;

	friend std::istream& operator>>(std::istream& in, AdjInfo& a)
	{
		std::string line;
		getline(in, line);

		// return if we've hit the end
		if(line == "")
			return in;
		
		StringVec fields = split(line, ',');
		assert(fields.size() == 4);

		VertexID from;
		VertexID to;
		int dir;
		bool comp;
	
		std::stringstream parser0(fields[0]);
		std::stringstream parser1(fields[1]);
		std::stringstream parser2(fields[2]);
		std::stringstream parser3(fields[3]);

		parser0 >> from;
		parser1 >> to;
		parser2 >> dir;
		parser3 >> comp;

		Edge tmp(from, to, (EdgeDir)dir, (EdgeComp)comp, a.overlap);
		a.e = tmp;
		return in;
	}
};


void loadContigVertices(SeqGraph& graph, int kmer, std::string filename)
{
	std::ifstream file(filename.c_str());
	assert(file.is_open());
	Contig c;
	(void)kmer;
	while(readCAF(file,c))
	{
		SeqVertex* pSV = new SeqVertex(c.getID(), c.getSequence());
		graph.addVertex(pSV);
	}
}

void loadContigEdges(int overlap, SeqGraph& graph, std::string filename)
{
	std::ifstream file(filename.c_str());
	assert(file.is_open());
	AdjInfo a(overlap);
	while(file >> a)
	{
		graph.addEdge(a.e);
	}
}

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string adjFile(argv[optind++]);
	
	SeqGraph sg;

	// Load verts and edges
	loadContigVertices(sg, opt::k, contigFile);
	loadContigEdges(opt::k - 1, sg, adjFile);

	//sg.stats();
	sg.validate();

	VertexIDVec nbVerts = sg.getNonBranchingVertices();
	for(VertexIDVec::iterator iter = nbVerts.begin(); iter != nbVerts.end(); ++iter)
	{
		std::cout << *iter << "\t" << "UNIQUE\n";
	}
	//sg.simplify();
	//sg.stats();
	//sg.validate();
}


void test()
{
	SeqGraph sg;
	Vertex* pV0 = new Vertex("0");
	Vertex* pV1 = new Vertex("1");
	Vertex* pV2 = new Vertex("2");
	Vertex* pV3 = new Vertex("3");
	Vertex* pV4 = new Vertex("4");
	sg.addVertex(pV0);
	sg.addVertex(pV1);
	sg.addVertex(pV2);
	sg.addVertex(pV3);
	sg.addVertex(pV4);

	sg.flip("4");
	sg.simplify();
	//sg.flip(4);

	sg.validate();
	sg.writeDot("blah");
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

