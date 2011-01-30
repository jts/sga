//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldAlgorithms - Collection of algorithms
// for computing linear scaffolds from a scaffold graph
//
#ifndef SCAFFOLDALGORITHMS_H
#define SCAFFOLDALGORITHMS_H

#include "ScaffoldSearch.h"
#include "ScaffoldWalk.h"
#include "ScaffoldGraph.h"

typedef std::vector<ScaffoldVertexPtrVector> ScaffoldConnectedComponents;

namespace ScaffoldAlgorithms
{
    // Construct scaffolds from the given graph
    void makeScaffolds(ScaffoldGraph* pGraph);

    // Compute the connected components of the scaffold graph
    void connectedComponents(ScaffoldGraph* pGraph, ScaffoldConnectedComponents& outComponents);

    // Compute the terminal vertices in the given connected component
    // A terminal vertex is one that has a connection in at most one direction
    // If the connected component forms a simple loop, this is possibly empty.
    // This case would signal an inconsistency in the scaffold graph 
    void computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector& component, 
                                               ScaffoldVertexPtrVector& terminals);

    // Compute the layout of a component of the scaffold graph starting from pVertex
    // to all the terminal vertices in the component.
    // Precondition: the connected component containing pVertex cannot have a cycle
    void computeLayout(ScaffoldVertex* pVertex, ScaffoldWalkVector& outWalks);

};

#endif
