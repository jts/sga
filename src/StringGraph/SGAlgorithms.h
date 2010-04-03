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
#include "SGUtil.h" 

namespace SGAlgorithms
{
	//
	// Overlap discovery algorithms
	//
	struct VertexOverlapPair
	{
		Vertex* pVertex;
		Overlap ovr;

		friend bool operator<(const VertexOverlapPair& a, const VertexOverlapPair& b)
		{
			return a.pVertex->getID() < b.pVertex->getID();
		}
	};

	typedef std::set<VertexOverlapPair> VertexOverlapSet;

	// Discover the complete set of overlaps for pVertex
	void findOverlapSet(const Vertex* pVertex, VertexOverlapSet& VOSet);

	// recursive function to discover overlaps of pX from pY
	void _discoverOverlaps(const Vertex* pX, const Vertex* pY, EdgeDir dir, const Overlap& ovrXY, 
	                      VertexOverlapSet& outSet);

	//
	// Overlap inference algorithms
	//
	// Infer an overlap from two edges
	// The input edges are between X->Y Y->Z
	// and the returned overlap is of the form X->Z
	Overlap inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);

	// Returns true if XZ has a non-zero length overlap
	bool hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);

	//
	// Construct an extended multioverlap for a vertex
	//
	MultiOverlap makeExtendedMultiOverlap(const Vertex* pVertex);

	//
	// Construct SeqTries from the extended overlap set
	//
	void makeExtendedSeqTries(const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie);
};

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
struct SGTransRedVisitor
{
	SGTransRedVisitor() {}
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

// Remodel the graph to infer missing edges or remove erroneous edges
struct SGRemodelVisitor
{
	SGRemodelVisitor() {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);
};

// 
struct SGErrorCorrectVisitor
{
	SGErrorCorrectVisitor() {}
	void previsit(StringGraph*) {}
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*) {}
};

// Infer missing edges in the graph
struct SGRealignVisitor
{
	struct Candidate
	{
		Candidate(Vertex* pv, const Overlap& o) : pEndpoint(pv), ovr(o) {}
		Vertex* pEndpoint;
		Overlap ovr;
	};
	typedef std::vector<Candidate> CandidateVector;

	SGRealignVisitor() {}
	void previsit(StringGraph* pGraph);
	bool visit(StringGraph* pGraph, Vertex* pVertex);
	void postvisit(StringGraph*);

	CandidateVector getMissingCandidates(StringGraph* pGraph, Vertex* pVertex, int minOverlap) const;
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
	void previsit(StringGraph*) {}
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
};

#endif
