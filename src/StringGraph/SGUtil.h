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
StringGraph* createSGFromOverlaps(std::string readFile, std::string overlapFile, std::string containFile);
void createVertices(StringGraph* pGraph, std::string readFile, const ContainMap& containments);
void createEdges(StringGraph* pGraph, std::string overlapFile, const ContainMap& containments);

//
std::string getOverhangString(const SeqCoord& sc, const std::string& seq);

#endif
