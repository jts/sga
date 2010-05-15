//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGAlgorithms.h"
#include "SGUtil.h"

#ifndef SGVISITORS_H
#define SGVISITORS_H

// Visit each node, writing it to a file as a fasta record
struct SGFastaVisitor
{
    // constructor
    SGFastaVisitor(std::string filename) : m_fileHandle(filename.c_str()) {}
    ~SGFastaVisitor() { m_fileHandle.close(); }

    // functions
    void previsit(StringGraph* /*pGraph*/) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* /*pGraph*/) {}

    // data
    std::ofstream m_fileHandle;
};

// Visit each node and write the overlaps to the specified file
struct SGOverlapWriterVisitor
{
    SGOverlapWriterVisitor(std::string filename) : m_fileHandle(filename.c_str()) {}
    ~SGOverlapWriterVisitor() { m_fileHandle.close(); }

    // functions
    void previsit(StringGraph* /*pGraph*/) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* /*pGraph*/) {}

    // data
    std::ofstream m_fileHandle;
};

// Run the Myers transitive reduction algorithm on each node
struct SGTransitiveReductionVisitor
{
    SGTransitiveReductionVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int marked_verts;
    int marked_edges;
};

// Remove contained vertices from the graph
struct SGContainRemoveVisitor
{
    SGContainRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
};

// Remove substring contained vertices from the graph
struct SGSubstringRemoveVisitor
{
    SGSubstringRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
};

// Validate that the graph does not contain
// any extra edges or missing irreducible edges
struct SGValidateStructureVisitor
{
    SGValidateStructureVisitor() {}
    void previsit(StringGraph*) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

// Remodel the graph to infer missing edges or remove erroneous edges
struct SGRemodelVisitor
{
    SGRemodelVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    double m_remodelER;
};

// 
struct SGErrorCorrectVisitor
{
    SGErrorCorrectVisitor() {}
    void previsit(StringGraph*) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

// Compute edge summary statistics 
struct SGEdgeStatsVisitor
{
    struct Candidate
    {
        Candidate(Vertex* pv, const Overlap& o) : pEndpoint(pv), ovr(o) {}
        Vertex* pEndpoint;
        Overlap ovr;
    };
    typedef std::vector<Candidate> CandidateVector;
    typedef std::map<int, int> IntIntMap;
    typedef std::map<int, IntIntMap> CountMatrix;

    //
    SGEdgeStatsVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    CandidateVector getMissingCandidates(StringGraph* pGraph, Vertex* pVertex, int minOverlap) const;
    void addOverlapToCount(int ol, int nd, CountMatrix& matrix);
    void printCounts(CountMatrix& matrix);

    //
    CountMatrix foundCounts;
    CountMatrix missingCounts;
    int maxDiff;
    int minOverlap;
    int maxOverlap;
};

// Detect whether vertices are dead ends and mark them for removal
struct SGTrimVisitor
{
    SGTrimVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int num_island;
    int num_terminal;
    int num_contig;
};

// Detect and remove duplicate edges
struct SGDuplicateVisitor
{
    SGDuplicateVisitor() {}
    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

// Detect small island vertices and removal them
struct SGIslandVisitor
{
    SGIslandVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);
};


// Detect whether vertices are bubbles and mark them for removal
struct SGBubbleVisitor
{
    SGBubbleVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);
    int num_bubbles;
};

// Compile summary statistics for the graph
struct SGGraphStatsVisitor
{
    SGGraphStatsVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int num_terminal;
    int num_island;
    int num_monobranch;
    int num_dibranch;
    int num_transitive;
    int num_edges;
    int num_vertex;
    size_t sum_edgeLen;
};

#endif
