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

// Compute the connected components of the graph
void ScaffoldAlgorithms::connectedComponents(ScaffoldGraph* pGraph)
{
    ScaffoldVertexPtrVector allVertices = pGraph->getAllVertices();

    ScaffoldConnectedComponents connectedComponents;
    ScaffoldSearchTree::connectedComponents(allVertices, connectedComponents);

    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        ScaffoldVertexPtrVector terminals;
        computeTerminalsForConnectedComponent(connectedComponents[i], terminals);
    }
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
// from a particular vertex.
// Precondition: The connected component containing pVertex does not contain
// a cycle.
void ScaffoldAlgorithms::computeLayout(ScaffoldVertex* pStartVertex)
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
        {
            walk.addEdge(*rIter);
        }
        std::cout << "walk to " << (*terminalIter)->getID() << " " << walk.getContigLengthSum() << " " << walk.getGapSum() << "\n";
        //walk.print();

    }
}

