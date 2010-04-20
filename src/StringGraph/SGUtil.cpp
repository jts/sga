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

StringGraph* SGUtil::loadASQG(const std::string& filename, const unsigned int minOverlap, 
                              bool allowContainments)
{
	// Initialize graph
	StringGraph* pGraph = new StringGraph;

	std::istream* pReader = createReader(filename);

	int stage = 0;
	int line = 0;
	std::string recordLine;
	while(getline(*pReader, recordLine))
	{
		ASQG::RecordType rt = ASQG::getRecordType(recordLine);
		switch(rt)
		{
			case ASQG::RT_HEADER:
			{
				if(stage != 0)
				{
					std::cerr << "Error: Unexpected header record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::HeaderRecord headerRecord(recordLine);
				(void)headerRecord; // do nothing with the header at this point
				break;
			}
			case ASQG::RT_VERTEX:
			{
				// progress the stage if we are done the header
				if(stage == 0)
					stage = 1;

				if(stage != 1)
				{
					std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::VertexRecord vertexRecord(recordLine);
				const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

				// If the substring tag is set and equals zero, the vertex is a substring
				// of some other vertex and should be skipped
				if(!ssTag.isInitialized() || ssTag.get() == 0)
				{
					pGraph->addVertex(new Vertex(vertexRecord.getID(), vertexRecord.getSeq()));
				}
				break;
			}
			case ASQG::RT_EDGE:
			{
				if(stage == 1)
					stage = 2;
				
				if(stage != 2)
				{
					std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
					exit(EXIT_FAILURE);
				}

				ASQG::EdgeRecord edgeRecord(recordLine);
				const Overlap& ovr = edgeRecord.getOverlap();
				if(!allowContainments && ovr.match.isContainment())
				{
					// Mark the contained read for subsequent removal
					const std::string& containedID = ovr.getContainedID();

					Vertex* pVertex = pGraph->getVertex(containedID);
					if(pVertex != NULL)
					{
						pVertex->setColor(GC_BLACK);
					}
				}
				else
				{
					// Add the edge to the graph
					if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
						SGUtil::createEdges(pGraph, ovr, allowContainments);
				}
				break;
			}
		}
		++line;
	}

	// Done loading the ASQG file, remove containment vertices if necessary and validate the graph
	pGraph->sweepVertices(GC_BLACK);

	// Remove any duplicate edges
	SGDuplicateVisitor dupVisit;
	pGraph->visit(dupVisit);

	delete pReader;
	return pGraph;
}

// add edges to the graph for the given overlap
Edge* SGUtil::createEdges(StringGraph* pGraph, const Overlap& o, bool allowContained)
{
	// Initialize data and perform checks
	Vertex* pVerts[2];
	EdgeComp comp = (o.match.isRC()) ? EC_REVERSE : EC_SAME;

	for(size_t idx = 0; idx < 2; ++idx)
	{
		pVerts[idx] = pGraph->getVertex(o.id[idx]);

		// If one of the vertices is not in the graph, skip this edge
		// This can occur if one of the verts is a strict substring of some other vertex so it will
		// never be added to the graph
		if(pVerts[idx] == NULL)
			return NULL;
		if(!allowContained)
			assert(!o.match.coord[idx].isContained() && o.match.coord[idx].isExtreme());
	}

	// Allocate the edges
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
