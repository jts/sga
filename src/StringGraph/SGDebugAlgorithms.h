//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
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
	int	num_trusted;
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

	void previsit(StringGraph*) {}
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*) {}

	// Data
	StringGraph* m_pCompareGraph;
};


#endif
