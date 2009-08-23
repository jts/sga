//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps
//
#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H

#include "Bigraph.h"
#include "Contig.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

//
// Implementations of data structures for the bigraph
//
struct SEData
{
	// constructor
	SEData(std::string s) : seq(s) {}

	// functions
	void flip();
	std::string getLabel() const;

	// data
	std::string seq;

};

typedef Edge StringEdge;

struct SVData
{
	// constructor
	SVData(std::string i, std::string s) : id(i), seq(s), readCount(1) {}

	// functions
	void merge(const SVData& other, const StringEdge& edge);

	// data
	std::string id;
	std::string seq;
	size_t readCount;
};

typedef Vertex StringVertex;
typedef Bigraph StringGraph;

// Visit functors
struct SGFastaVisitor
{
	// constructor
	SGFastaVisitor(std::string filename) : m_fileHandle(filename.c_str()) {}
	~SGFastaVisitor() { m_fileHandle.close(); }

	// functions
	bool visit(StringGraph* pGraph, StringVertex* pVertex);

	// data
	std::ofstream m_fileHandle;
};

// functions
Sequence flip(const Sequence& s);
StringGraph* createStringGraph(std::string readFile, std::string overlapFile);
std::string getOverhangString(const SeqCoord& sc, const std::string& seq);

#endif
