//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGPairedAlgorithms - Collection of algorithms
// using paired end data in string graphs
//
#ifndef SGPAIREDALGORITHMS_H
#define SGPAIREDALGORITHMS_H

#include "Bigraph.h"
#include "SGUtil.h" 

// Connect the paired reads in the string graph
struct SGVertexPairingVisitor
{
    SGVertexPairingVisitor() {}
    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int num_paired;
    int num_unpaired;
};

// Build the paired end trust network
struct SGPETrustVisitor
{
    SGPETrustVisitor() {}
    void previsit(StringGraph*) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

// Remove untrusted, conflicted edges from the graph
// Experimental.
struct SGPEConflictRemover
{
    SGPEConflictRemover() {}
    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int num_same;
    int num_diff;

};

// Visit each node and output the overlap between each linked edge and their pairs
struct SGPairedOverlapVisitor
{
    SGPairedOverlapVisitor() {}
    void previsit(StringGraph*) {} 
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

#endif
