//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#include "SGUtil.h"
#include "SeqReader.h"
#include "SGAlgorithms.h"

ContainMap::ContainMap(std::string file)
{
	std::ifstream reader(file.c_str());
	checkFileHandle(reader, file);
	Overlap o;
	while(reader >> o)
	{
		size_t idx = o.getContainedIdx();
		add(o.id[idx], o.id[1 - idx]);
	}
	reader.close();
}

void ContainMap::add(std::string s1, std::string s2)
{
	m_data.insert(std::make_pair(s1, s2));
}

bool ContainMap::isContained(std::string s) const
{
	StrStrMap::const_iterator iter = m_data.find(s);
	return iter != m_data.end();
}

//
std::string ContainMap::getContainer(std::string s) const
{
	StrStrMap::const_iterator iter = m_data.find(s);
	if(iter != m_data.end())
		return iter->second;
	else
		return "";
}

// Construct a string graph from overlaps
StringGraph* loadStringGraph(std::string readFile, std::string overlapFile, std::string containFile, const unsigned int minOverlap, bool allowContainments)
{
	// Initialize the string graph
	StringGraph* pGraph = new StringGraph;

	// Load the containment mappings
	ContainMap containments(containFile);
	
	// Create the graph
	loadVertices(pGraph, readFile, containments, allowContainments);
	loadEdges(pGraph, overlapFile, containments, minOverlap, allowContainments);

	if(allowContainments)
		loadEdges(pGraph, containFile, containments, minOverlap, allowContainments);

	// Validate the graph and ensure that there are no duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	return pGraph;
}

//
void loadVertices(StringGraph* pGraph, std::string readFile, const ContainMap& containments, bool allowContainments)
{
	(void)containments;
	// Add the reads as the vertices
	SeqReader reader(readFile);
	SeqRecord sr;
	while(reader.get(sr))
	{
		if(allowContainments || !containments.isContained(sr.id))
		{
			pGraph->addVertex(new Vertex(sr.id, sr.seq.toString()));
		}
	}
}

//
void loadEdges(StringGraph* pGraph, std::string overlapFile, const ContainMap& containments, 
               const unsigned int minOverlap, bool allowContainments)
{
	// Add the overlaps as edges
	std::ifstream overlapReader(overlapFile.c_str());
	checkFileHandle(overlapReader, overlapFile);

	Overlap o;
	while(overlapReader >> o)
	{
	    if((!allowContainments && (containments.isContained(o.id[0]) || containments.isContained(o.id[1])))
			|| o.match.getMaxOverlapLength() < (int)minOverlap)
		{
			continue;
		}
		createEdges(pGraph, o, allowContainments);
	}
	overlapReader.close();
}

// add edges to the graph for the given overlap
Edge* createEdges(StringGraph* pGraph, const Overlap& o, bool allowContained)
{
	// Initialize data and perform checks
	Vertex* pVerts[2];
	EdgeComp comp = (o.match.isRC()) ? EC_REVERSE : EC_SAME;

	for(size_t idx = 0; idx < 2; ++idx)
	{
		pVerts[idx] = pGraph->getVertex(o.id[idx]);
		assert(pVerts[idx]);		
		if(!allowContained)
			assert(!o.match.coord[idx].isContained() && o.match.coord[idx].isExtreme());
	}

	// Allocated the edges
	Edge* pEdges[2];
	for(size_t idx = 0; idx < 2; ++idx)
	{
		EdgeDir dir = o.match.coord[idx].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE;
		const SeqCoord& coord = o.match.coord[idx];
		pEdges[idx] = new Edge(pVerts[1 - idx], dir, comp, coord);
	}

	pEdges[0]->setTwin(pEdges[1]);
	pEdges[1]->setTwin(pEdges[0]);
	
	bool isContainment = o.match.isContainment();
	pGraph->addEdge(pVerts[0], pEdges[0]);
	pGraph->addEdge(pVerts[1], pEdges[1]);

	if(isContainment)
		pGraph->setContainmentFlag(true);

	return pEdges[0];
}
