//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#ifndef SGUTIL_H
#define SGUTIL_H

#include "StringGraph.h"

// typedefs
typedef std::map<std::string, std::string> StrStrMap;

// data
struct ContainMap
{
	public:
		// Constructors
		ContainMap(std::string file);

		// functions
		void add(std::string s1, std::string s2);
		bool isContained(std::string s) const;
		std::string getContainer(std::string s) const;

	private:
		StrStrMap m_data;
};

// functions

// string graph creation
StringGraph* loadStringGraph(std::string readFile, std::string overlapFile, std::string containFile, const unsigned int minOverlap);
void loadVertices(StringGraph* pGraph, std::string readFile, const ContainMap& containments);
void loadEdges(StringGraph* pGraph, std::string overlapFile, const ContainMap& containments, const unsigned int minOverlap);

// Create the edges described by the overlap. Returns a pointer to the edge
// from the first entry of the overlap to the second. 
StringEdge* createEdges(StringGraph* pGraph, const Overlap& o);

// Return the index (0 or 1) of the CONTAINED vertex (the discarded vertex of a containment)
size_t getContainedIdx(const Overlap& o);


//
std::string getOverhangString(const SeqCoord& sc, const std::string& seq);

#endif
