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
void ScaffoldAlgorithms::computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector component, 
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

