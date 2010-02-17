//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGDebugAlgorithms - Methods used for the 
// development of algorithms.
//
#include "SGDebugAlgorithms.h"

//
// SGDebugEdgeClassificationVisitor - Collect statistics about the graph
// using debug information about simulated reads
//
void SGDebugEdgeClassificationVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_good = 0;
	num_bad = 0;
	num_conflicted = 0;
	num_trusted = 0;
	num_nottrusted = 0;
}

//
bool SGDebugEdgeClassificationVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	assert(false && "Todo: reimplement");
	int curr_pos = 0; //pVertex->dbg_position;
	EdgePtrVec edges = pVertex->getEdges();

	if(pVertex->countEdges(ED_SENSE) > 1 || pVertex->countEdges(ED_ANTISENSE) > 1)
		num_conflicted++;

	for(size_t i = 0; i < edges.size(); ++i)
	{
		int edge_pos = 0; //edges[i]->getEnd()->dbg_position;
		int distance = abs(curr_pos - edge_pos);
		int overlap_len = edges[i]->getMatchLength();
		int sum = distance + overlap_len;
		
		if(sum == (int)pVertex->getSeq().size())
			++num_good;
		else
			++num_bad;

		if(edges[i]->isTrusted)
			num_trusted++;
		else
			num_nottrusted++;
	}
	return false;
}

// Remove all the marked edges
void SGDebugEdgeClassificationVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("Num good: %d Num bad: %d Num conflicted: %d Num trusted: %d Num not trusted: %d\n", 
	        num_good, num_bad, num_conflicted, num_trusted, num_nottrusted);
}

//
// SGDebugGraphCompareVisitor
// Compare the visited graph to the graph loaded
// in the constructor
//
SGDebugGraphCompareVisitor::SGDebugGraphCompareVisitor(std::string readsFile)
{
	std::string prefix = stripFilename(readsFile);
	std::cout << "Loading graph of " << readsFile << " for comparison\n";
	m_pCompareGraph = loadStringGraph(readsFile, prefix + ".ovr", prefix + ".ctn", 0);
}

//
SGDebugGraphCompareVisitor::~SGDebugGraphCompareVisitor()
{
	delete m_pCompareGraph;
	m_pCompareGraph = 0;
}

//
bool SGDebugGraphCompareVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	// Retreive the vertex in the comparison graph
	Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
	if(pCompareVertex == NULL)
	{
		std::cerr << "Compare vertex " << pVertex->getID() << " not found!\n";
		return false;
	}
	EdgePtrVec compareEdges = pCompareVertex->getEdges();

	for(size_t i = 0; i < compareEdges.size(); ++i)
	{
		Edge* pCompareEdge = compareEdges[i];
		EdgeDesc ed = pCompareEdge->getDesc();
		bool is_missing = !pVertex->hasEdge(ed);

		//if(!pVertex->hasEdge(ed))
		{
			//std::cout << "Graph is missing edge " << ed << "\n";

			// Get the expected match between the sequences
			Match exp_match = pCompareEdge->getMatch();

			// Get the sequences of the endpoints for the current graph
			Vertex* pVertex2 = pGraph->getVertex(pCompareEdge->getEndID());

			std::string match1 = exp_match.coord[0].getSubstring(pVertex->getSeq());
			std::string match2 = exp_match.coord[1].getSubstring(pVertex2->getSeq());
			if(exp_match.isRC())
				match2 = reverseComplement(match2);
			
			std::string ms = is_missing ? "MISSING" : "FOUND";
			int num_diffs = countDifferences(match1, match2, match1.size());
			double mm_rate = double(num_diffs) / double(match1.size());
			(void)mm_rate;
			printf("%s\t%d\t%zu\t%lf\n", ms.c_str(), num_diffs, match1.size(), mm_rate);
			//std::cout << "SEQ1: " << pVertex->getSeq() << "\n";
			//std::cout << "SEQ2: " << pVertex2->getSeq() << "\n";
			//std::cout << "MATCH1: " << match1 << "\n";
			//std::cout << "MATCH2: " << match2 << "\n";
		}
	}

	return false;
}
