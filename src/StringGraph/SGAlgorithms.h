//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#ifndef SGALGORITHMS_H
#define SGALGORITHMS_H

#include "Bigraph.h"
#include "SGUtil.h"
#include <queue>

namespace SGAlgorithms
{
	//
	// Overlap discovery algorithms
	//

	// Structure used for iterative exploration of the graph
	// The exploration starts at some vertex X, each element
	// holds a possible overlap between X and some other vertex Y
	struct ExploreElement
	{
		ExploreElement(const EdgeDesc& e, const Overlap& o) : ed(e), ovr(o) {}
		EdgeDesc ed;
		Overlap ovr;
	};

	// Comparison operator used to compare ExploreElements
	// by the length of the overlap on vertex X
	struct CompareExploreElemOverlapLength
	{
		bool operator()(const ExploreElement& elemXY, const ExploreElement& elemXZ)
		{
			return elemXY.ovr.match.coord[0].length() < elemXZ.ovr.match.coord[0].length();
		}
	};

	typedef std::queue<ExploreElement> ExploreQueue;
	typedef std::priority_queue<ExploreElement, 
	                            std::vector<ExploreElement>, 
								CompareExploreElemOverlapLength> ExplorePriorityQueue;

	//
	typedef std::map<EdgeDesc, Overlap> EdgeDescOverlapMap;
	typedef std::set<EdgeDesc> EdgeDescSet;

	// Remodel the edges of pVertex by finding any new irreducible edges
	// that may need to be added if the vertex at the end of pEdge is removed
	// from the graph
	void remodelVertexAfterRemoval(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

	// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
	void remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

	// Add the neighbors of the endpoint of edXY to the explore queue if they overlap pX
	void enqueueEdges(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
                      ExplorePriorityQueue& outQueue, EdgeDescSet& seenEdges, 
					  EdgeDescOverlapMap* pExclusionSet);

	// Add overlaps to pX inferred from the edges of pY to outMap
	void addOverlapsToSet(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, EdgeDescOverlapMap& outMap);

	// Discover the complete set of overlaps for pVertex
	void findOverlapMap(const Vertex* pVertex, EdgeDescOverlapMap& outMap);

	// recursive function to discover overlaps of pX from pY
	void _discoverOverlaps(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
	                       EdgeDescOverlapMap& outMap);

	// Calculate the error rate between the two vertices
	double calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY);

	//
	// Overlap inference algorithms
	//
	// Infer an overlap from two edges
	// The input overlaps are between X->Y Y->Z
	// and the returned overlap is X->Z
	Overlap inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);
	EdgeDesc inferTransitiveEdgeDesc(const EdgeDesc& edXY, const EdgeDesc& edYZ);

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
