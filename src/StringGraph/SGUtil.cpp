//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#include "SGUtil.h"
#include "SeqReader.h"

ContainMap::ContainMap(std::string file)
{
	std::ifstream reader(file.c_str());
	checkFileHandle(reader, file);
	Overlap o;
	while(reader >> o)
	{
		size_t idx = getContainedIdx(o);
		add(o.id[idx], o.id[1 - idx]);
	}
	reader.close();
}

void ContainMap::add(std::string s1, std::string s2)
{
	/*
	if(isContained(s1))
		std::cerr << "Warning: " << s1 << " is contained in: " << getContainer(s1) << "\n";
	if(isContained(s2))
		std::cerr << "Warning: " << s2 << " is contained in: " << getContainer(s2) << "\n";
	*/
	m_data.insert(std::make_pair(s1, s2));
}

bool ContainMap::isContained(std::string s) const
{
	StrStrMap::const_iterator iter = m_data.find(s);
	return iter != m_data.end();
}

std::string ContainMap::getContainer(std::string s) const
{
	StrStrMap::const_iterator iter = m_data.find(s);
	if(iter != m_data.end())
		return iter->second;
	else
		return "";
}

// Return the index of the CONTAINED vertex
size_t getContainedIdx(const Overlap& o)
{
	// The verts are mutually contained, return the lexographically lower id
	if(o.match.coord[0].isContained() && o.match.coord[1].isContained())
	{
		if(o.id[0] < o.id[1])
			return 1;
		else
			return 0;
	}
	else if(o.match.coord[0].isContained())
	{
		return 0;
	}
	else
	{
		assert(o.match.coord[1].isContained());
		return 1;
	}
}

// Construct a string graph from overlaps
StringGraph* loadStringGraph(std::string readFile, std::string overlapFile, std::string containFile, const unsigned int minOverlap)
{
	// Initialize the string graph
	StringGraph* pGraph = new StringGraph;

	// Load the containment mappings
	ContainMap containments(containFile);
	
	// Create the graph
	loadVertices(pGraph, readFile, containments);
	loadEdges(pGraph, overlapFile, containments, minOverlap);
	return pGraph;
}

void loadVertices(StringGraph* pGraph, std::string readFile, const ContainMap& containments)
{
	// Add the reads as the vertices
	SeqReader reader(readFile);
	SeqItem si;
	while(reader.get(si))
	{
		if(!containments.isContained(si.id))
		{
			pGraph->addVertex(new StringVertex(si.id, si.seq.toString()));
		}
	}
}

void loadEdges(StringGraph* pGraph, std::string overlapFile, const ContainMap& containments, const unsigned int minOverlap)
{
	// Add the overlaps as edges
	std::ifstream overlapReader(overlapFile.c_str());
	checkFileHandle(overlapReader, overlapFile);

	Overlap o;
	while(overlapReader >> o)
	{
		if(containments.isContained(o.id[0]) || containments.isContained(o.id[1]) || o.match.getMaxOverlapLength() < (int)minOverlap)
			continue;
		createEdges(pGraph, o);
	}
	overlapReader.close();
}

// add edges to the graph for the given overlap
StringEdge* createEdges(StringGraph* pGraph, const Overlap& o)
{
	// Initialize data and perform checks
	StringVertex* pVerts[2];
	EdgeComp comp = (o.match.isRC()) ? EC_REVERSE : EC_SAME;

	for(size_t idx = 0; idx < 2; ++idx)
	{
		pVerts[idx] = static_cast<StringVertex*>(pGraph->getVertex(o.id[idx]));
		assert(pVerts[idx]);		
		// Ensure the reads are not identical
		assert(!o.match.coord[idx].isContained() && o.match.coord[idx].isExtreme());
	}

	// Allocated the edges
	StringEdge* pEdges[2];
	for(size_t idx = 0; idx < 2; ++idx)
	{
		EdgeDir dir = o.match.coord[idx].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE;
		const SeqCoord& coord = o.match.coord[idx];
		pEdges[idx] = new StringEdge(pVerts[idx], pVerts[1 - idx], dir, comp, coord, o.match.getNumDiffs());
	}

	pEdges[0]->setTwin(pEdges[1]);
	pEdges[1]->setTwin(pEdges[0]);

	pGraph->addEdge(pEdges[0]);
	pGraph->addEdge(pEdges[1]);

	return pEdges[0];
}
