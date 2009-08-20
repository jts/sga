//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps
//
#include "StringGraph.h"
#include "SeqReader.h"

//
StringGraph* createStringGraph(std::string readFile, std::string overlapFile)
{
	StringGraph* pGraph = new StringGraph;

	// Add the reads as the vertices
	SeqReader reader(readFile);
	SeqItem si;
	while(reader.get(si))
		pGraph->addVertex(new StringVertex(si.id, si));

	// Add the overlaps as edges
	std::ifstream overlapReader(overlapFile.c_str());
	Overlap o;
	while(overlapReader >> o)
	{
		std::cout << "Read: " << o << "\n";
	}
	overlapReader.close();
	return pGraph;
}

