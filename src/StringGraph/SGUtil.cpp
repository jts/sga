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
                const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
                if(overlapTag.isInitialized())
                    pGraph->setMinOverlap(overlapTag.get());


                const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
                if(errorRateTag.isInitialized())
                    pGraph->setErrorRate(errorRateTag.get());

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

    bool isContainment = o.match.isContainment();
    assert(allowContained || !isContainment);

    for(size_t idx = 0; idx < 2; ++idx)
    {
        pVerts[idx] = pGraph->getVertex(o.id[idx]);

        // If one of the vertices is not in the graph, skip this edge
        // This can occur if one of the verts is a strict substring of some other vertex so it will
        // never be added to the graph
        if(pVerts[idx] == NULL)
            return NULL;
        assert(o.match.coord[idx].isExtreme());
    }
    if(!isContainment)
    {
        Edge* pEdges[2];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            EdgeDir dir = o.match.coord[idx].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE;
            const SeqCoord& coord = o.match.coord[idx];
            pEdges[idx] = new Edge(pVerts[1 - idx], dir, comp, coord);
        }

        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);
        
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[1], pEdges[1]);
        return pEdges[0];
    }
    else
    {
        // Contained edges don't have a direction, they can be travelled from
        // one vertex to the other in either direction. Hence, we 
        // add two edges per vertex. Later during the contain removal
        // algorithm this is important to determine transitivity
        Edge* pEdges[4];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            const SeqCoord& coord = o.match.coord[idx];
            pEdges[idx] = new Edge(pVerts[1 - idx], ED_SENSE, comp, coord);
            pEdges[idx + 2] = new Edge(pVerts[1 - idx], ED_ANTISENSE, comp, coord);
        }
        
        // Twin the edges and add them to the graph
        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);

        pEdges[2]->setTwin(pEdges[3]);
        pEdges[3]->setTwin(pEdges[2]);
    
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[0], pEdges[2]);

        pGraph->addEdge(pVerts[1], pEdges[1]);
        pGraph->addEdge(pVerts[1], pEdges[3]);

        pGraph->setContainmentFlag(true);
        return pEdges[0];
    }
}
