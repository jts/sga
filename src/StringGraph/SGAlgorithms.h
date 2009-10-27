//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#ifndef SGALGORITHMS_H
#define SGALGORITHMS_H

#include "Bigraph.h"
#include "StringGraph.h"

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

// Run the Myers transitive reduction algorithm on each node
struct SGTransRedVisitor
{
	SGTransRedVisitor() {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);

	int marked_verts;
	int marked_edges;
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

// Detect whether vertices are bubbles and mark them for removal
struct SGBubbleVisitor
{
	SGBubbleVisitor() {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);
	int num_bubbles;
};

// Discover variants through a local realignment
struct SGVariantVisitor
{
	SGVariantVisitor() : m_fileHandle("inferred.ovr") {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);

	std::ofstream m_fileHandle;						 
};

// Infer an error rate for each node are mark nodes about the error threshold for removal
struct SGErrorRemovalVisitor
{
	SGErrorRemovalVisitor(double er) : m_maxErrorRate(er) {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);
	double m_maxErrorRate;
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
};

#endif
