//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGDebugAlgorithms - Methods used for the 
// development of algorithms.
//
#ifndef SGDEBUGALGORITHMS_H
#define SGDEBUGALGORITHMS_H

#include "Bigraph.h"
#include "SGUtil.h"

// Classify edges from simulated data
// as good (correct) or bad (incorrect) based 
// on positions encoded in the read name.
struct SGDebugEdgeClassificationVisitor
{
    SGDebugEdgeClassificationVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int getNumGood() { return num_good; }
    int getNumBad() { return num_bad; }

    // Data
    int num_good;
    int num_bad;
    int num_conflicted;
    int    num_trusted;
    int num_nottrusted;
};

//
// Compare the visited graph to the graph loaded
// in the constructor
//
struct SGDebugGraphCompareVisitor
{
    SGDebugGraphCompareVisitor(std::string readsFile);
    ~SGDebugGraphCompareVisitor();

    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    void summarize(StringGraph* pGraph, Vertex* pVertex);
    void showMissing(StringGraph* pGraph, Vertex* pVertex);

    void compareErrorRates(StringGraph* pGraph, Vertex* pVertex);
    void compareTransitiveGroups(StringGraph* pGraph, Vertex* pVertex);
    
    // Data
    StringGraph* m_pCompareGraph;
    int m_numContained;
    int m_numFound;
    int m_numMissingNull;
    int m_numMissing;
    int m_numWrong;
    int m_numClosed;
    int m_numOpen;

    bool m_showMissing;
};


#endif
