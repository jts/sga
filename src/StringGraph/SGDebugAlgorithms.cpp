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

bool SGDebugEdgeClassificationVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	int curr_pos = pVertex->dbg_position;
	EdgePtrVec edges = pVertex->getEdges();

	if(pVertex->countEdges(ED_SENSE) > 1 || pVertex->countEdges(ED_ANTISENSE) > 1)
		num_conflicted++;

	for(size_t i = 0; i < edges.size(); ++i)
	{
		int edge_pos = edges[i]->getEnd()->dbg_position;
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
	printf("Num good: %d Num bad: %d Num conflicted: %d Num trusted: %d Num not trusted: %d\n", num_good, num_bad, num_conflicted, num_trusted, num_nottrusted);
}
