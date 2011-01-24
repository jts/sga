//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSearch - Find walks through a scaffold graph
//
#include "ScaffoldSearch.h"
#include "ScaffoldVertex.h"

struct ScaffoldDistanceFunction
{
    int operator()(const ScaffoldEdge* pEdge) const
    {
        return pEdge->getDistance();
    }
};
typedef GraphSearchTree<ScaffoldVertex, ScaffoldEdge, ScaffoldDistanceFunction> ScaffoldSearchTree;

void ScaffoldSearch::findVariantWalks(ScaffoldVertex* pX, 
                                      EdgeDir initialDir, 
                                      int maxDistance,
                                      size_t maxWalks, 
                                      ScaffoldWalkVector& outWalks)
{
    (void)maxWalks;
    ScaffoldSearchTree searchTree(pX, NULL, initialDir, maxDistance, 10000);

    // Iteravively perform the BFS using the search tree. After each step
    // we check if the search has collapsed to a single vertex.
    bool done = false;
    while(!done)
    {
        done = !searchTree.stepOnce();
        if(searchTree.wasSearchAborted())
            break;

        ScaffoldVertex* pCollapsedVertex;
        bool isCollapsed = searchTree.hasSearchConverged(pCollapsedVertex);
        if(isCollapsed)
        {
            assert(pCollapsedVertex != NULL);
            // pCollapsedVertex is common between all walks.
            // Check that the extension distance for any walk is no longer than maxDistance.
            // If this is the case, we truncate all walks so they end at pCollapsedVertex and return 
            // all the non-redundant walks in outWalks
            
            VertexID iLastID = pCollapsedVertex->getID();
            std::cout << "Walk collapsed at " << iLastID << "\n";

            // initial walks
            SEPVVector edgeWalks;
            searchTree.buildWalksToAllLeaves(edgeWalks);

            ScaffoldWalkVector walkVector;
            convertEdgeVectorsToScaffoldWalk(edgeWalks, walkVector); 
            return;           
        }
    }

    // no collapsed walk found, return empty set
    outWalks.clear();
}

void ScaffoldSearch::convertEdgeVectorsToScaffoldWalk(const SEPVVector& edgeWalks, 
                                                ScaffoldWalkVector& outWalks)
{
    // Convert the walks described by edges into SGWalks
    for(SEPVVector::const_iterator iter = edgeWalks.begin();
                                              iter != edgeWalks.end();
                                              ++iter)
    {
        outWalks.push_back(ScaffoldWalk(*iter));
    }
}
