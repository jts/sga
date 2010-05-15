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

namespace SGPairedAlgorithms
{
    // Find paths between the two vertices that are no longer than maxDistance
    void searchPaths(const Vertex* pX, const Vertex* pY, int maxDistance, PathVector& outPaths);
    std::string pathToString(const Vertex* pX, const Path& path);

    //
    EdgeDir getDirectionToPair(const std::string& id);

};

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

// Complete the sequence between a vertex and its pair
struct SGPairedPathResolveVisitor
{
    public:
        SGPairedPathResolveVisitor();
        ~SGPairedPathResolveVisitor();

        void previsit(StringGraph*);
        bool visit(StringGraph* pGraph, Vertex* pVertex);
        void postvisit(StringGraph*);

    private:
        std::ostream* m_pWriter;
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
