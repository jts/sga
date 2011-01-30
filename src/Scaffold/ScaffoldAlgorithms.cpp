//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldAlgorithms - Collection of algorithms
// for computing linear scaffolds from a scaffold graph
//
#include "ScaffoldAlgorithms.h"

// Build the scaffolds.
// The algorithm first finds the connected components of the graph then
// all the terminal vertices of each component. It then searches the
// components starting from the terminals for a walk that maximizes
// the size of the contigs in the primary scaffold.
void ScaffoldAlgorithms::makeScaffolds(ScaffoldGraph* pGraph)
{
    // Set every edge to be colored black
    // Edges that are kept in the scaffold will be colored white
    pGraph->setEdgeColors(GC_BLACK);

    ScaffoldConnectedComponents connectedComponents;
    ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);
 
    // Iterate over each connected component and linearize the scaffolds   
    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        WARN_ONCE("CHECK FOR CYCLES IN SCAFFOLD GRAPH");


        // Discover the terminal vertices in the component
        ScaffoldVertexPtrVector& component = connectedComponents[i];

        if(component.size() == 1) 
            continue; //nothing to do

        //
        ScaffoldVertexPtrVector terminalVertices;
        computeTerminalsForConnectedComponent(component, terminalVertices);

        // Construct a scaffold layout for each terminal vertex of the component
        // We select the layout that contains the greatest amount of sequence as
        // the initial layout of the scaffold
        size_t bestLayoutBases = 0;
        ScaffoldWalk bestWalk(NULL);

        for(size_t j = 0; j < terminalVertices.size(); j++)
        {
            ScaffoldWalkVector layouts;
            computeLayout(terminalVertices[j], layouts);

            for(size_t k = 0; k < layouts.size(); ++k)
            {
                size_t layoutSum = layouts[k].getContigLengthSum();

                if(layoutSum > bestLayoutBases)
                {
                    bestLayoutBases = layoutSum;
                    bestWalk = layouts[k];
                }
            }
        }

        // Ensure a valid layout is chosen
        assert(bestLayoutBases > 0 && bestWalk.getStartVertex() != NULL);

        // Remove every edge that is not a part of the chosen walk
        // Since every edge in the graph was initially colored black,
        // we color the kept edges white then remove every black edge
        // in a single pass later
        ScaffoldEdgePtrVector keptEdges = bestWalk.getEdges();
        for(size_t j = 0; j < keptEdges.size(); ++j)
        {
            keptEdges[j]->setColor(GC_WHITE);
            keptEdges[j]->getTwin()->setColor(GC_WHITE);
        }
    }

    // Remove all edges that aren't on a chosen layouts
    pGraph->deleteEdgesByColor(GC_BLACK);
}

// Compute the connected components of the graph
void ScaffoldAlgorithms::connectedComponents(ScaffoldGraph* pGraph, ScaffoldConnectedComponents& outComponents)
{
    ScaffoldVertexPtrVector allVertices = pGraph->getAllVertices();
    ScaffoldSearchTree::connectedComponents(allVertices, outComponents);
}

// Compute the terminal vertices of a connected component
void ScaffoldAlgorithms::computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector& component, 
                                                               ScaffoldVertexPtrVector& terminals)
{
    for(size_t i = 0; i < component.size(); ++i)
    {
        ScaffoldVertex* pVertex = component[i];
        size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
        size_t sCount = pVertex->getEdges(ED_SENSE).size();

        if(asCount == 0 || sCount == 0)
            terminals.push_back(pVertex);
    }

    std::cout << "CC: ";
    for(size_t i = 0; i < component.size(); ++i)
        std::cout << component[i]->getID() << " ";
    std::cout << "\n";

    std::cout << "Terminals: ";
    for(size_t i = 0; i < terminals.size(); ++i)
        std::cout << terminals[i]->getID() << " ";
    std::cout << "\n";
}

struct LayoutNode
{
    ScaffoldEdge* pEdge; // Current edge being extended from
    int distance; // Sum of the gap lengths to the endpoint of the edge
};

typedef std::map<ScaffoldVertex*, int> LayoutDistanceMap;
typedef std::map<ScaffoldVertex*, ScaffoldEdge*> LayoutEdgeMap;
typedef std::queue<LayoutNode> LayoutQueue;
typedef std::set<ScaffoldVertex*> LayoutTerminalSet;

// Compute a layout of a component of the scaffold graph starting
// from a particular vertex. Writes a path each reachable terminal vertex
// into outWalks.
// Precondition: The connected component containing pVertex does not contain
// a cycle.
void ScaffoldAlgorithms::computeLayout(ScaffoldVertex* pStartVertex, ScaffoldWalkVector& outWalks)
{
    // Detect the direction to begin constructing the component
    size_t asCount = pStartVertex->getEdges(ED_ANTISENSE).size();
    size_t sCount = pStartVertex->getEdges(ED_SENSE).size();
    assert(asCount == 0 || sCount == 0);
    if(asCount == 0 && sCount == 0)
        return; // nothing to do
    EdgeDir dir = (asCount > 0) ? ED_ANTISENSE : ED_SENSE;

    LayoutDistanceMap distanceMap;
    LayoutEdgeMap edgeMap;
    LayoutTerminalSet terminalSet;

    distanceMap[pStartVertex] = 0;
    edgeMap[pStartVertex] = NULL;

    // Enqueue the edges of this node to the search structure
    LayoutQueue visitQueue;
    ScaffoldEdgePtrVector startEdges = pStartVertex->getEdges(dir);
    for(size_t i = 0; i < startEdges.size(); ++i)
    {
        ScaffoldVertex* pZ = startEdges[i]->getEnd();
        LayoutNode initial;
        initial.pEdge = startEdges[i];
        initial.distance = startEdges[i]->getDistance();

        distanceMap[pZ] = initial.distance;
        edgeMap[pZ] = startEdges[i];

        visitQueue.push(initial);
    }

    while(!visitQueue.empty())
    {
        LayoutNode node = visitQueue.front();
        visitQueue.pop();
        ScaffoldEdge* pXY = node.pEdge;
        ScaffoldVertex* pY = pXY->getEnd();

        // Get the edges of the endpoint and add them to the queue
        // to explore they are closer than the current estimate
        EdgeDir yDir = !pXY->getTwin()->getDir();
        ScaffoldEdgePtrVector yEdges = pY->getEdges(yDir);
        if(yEdges.empty())
        {
            terminalSet.insert(pY);
        }

        // Add the edges of y to the queue if a new or shorter path has been found
        for(size_t i = 0; i < yEdges.size(); ++i)
        {
            ScaffoldEdge* pYZ = yEdges[i];
            ScaffoldVertex* pZ = pYZ->getEnd();

            int zDistance = node.distance + pYZ->getDistance();

            // Check if the gap distance to z is less than the current
            // If so, update the edge to z to be the current edge and re-queue
            // the edge
            LayoutDistanceMap::iterator distIter = distanceMap.find(pZ);
            if(distIter == distanceMap.end() || distIter->second > zDistance)
            {
                distanceMap[pZ] = zDistance;
                edgeMap[pZ] = pYZ;
                LayoutNode updateNode;
                updateNode.pEdge = pYZ;
                updateNode.distance = zDistance;
                visitQueue.push(updateNode);
            }
        }
    }

    // Construct paths from the starting node to the terminal vertices
    LayoutTerminalSet::iterator terminalIter = terminalSet.begin();
    for(; terminalIter != terminalSet.end(); ++terminalIter)
    {
        ScaffoldEdgePtrVector reverseEdges;
        ScaffoldVertex* pCurrent = *terminalIter;
        while(pCurrent != pStartVertex)
        {
            ScaffoldEdge* pEdge = edgeMap[pCurrent];
            reverseEdges.push_back(pEdge);
            assert(pEdge != NULL);
            pCurrent = pEdge->getStart();
        }

        // Build the real walk
        ScaffoldWalk walk(pStartVertex);
        ScaffoldEdgePtrVector::reverse_iterator rIter = reverseEdges.rbegin();
        for(; rIter != reverseEdges.rend(); ++rIter)
            walk.addEdge(*rIter);
        outWalks.push_back(walk);
    }
}

