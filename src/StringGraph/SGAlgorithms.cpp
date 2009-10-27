//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#include "SGAlgorithms.h"

// Visitor which outputs the graph in fasta format
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	StringVertex* pSV = static_cast<StringVertex*>(pVertex);
	m_fileHandle << ">" << pSV->getID() << " " <<  pSV->getSeq().length() << " " << pSV->getReadCount() << "\n";
	m_fileHandle << pSV->getSeq() << "\n";
	return false;
}

void SGTransRedVisitor::previsit(StringGraph* pGraph)
{
	// Set all the vertices in the graph to "vacant"
	pGraph->setColors(GC_WHITE);
	std::cout << "Running TR algorithm in EXACT mode\n";
	pGraph->sortVertexAdjLists();

	marked_verts = 0;
	marked_edges = 0;
}

// Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
bool SGTransRedVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	(void)pGraph;
	static const size_t FUZZ = 10; // see myers...

	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		//std::cout << pVertex->getID() << "," << dir << "\n";
		EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted

		if(edges.size() == 0)
			continue;

		for(size_t i = 0; i < edges.size(); ++i)
			(edges[i])->getEnd()->setColor(GC_GRAY);

		StringEdge* pLongestEdge = static_cast<StringEdge*>(edges.back());
		size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;
		
		// Stage 1
		for(size_t i = 0; i < edges.size(); ++i)
		{
			StringEdge* pVWEdge = static_cast<StringEdge*>(edges[i]);
			StringVertex* pWVert = static_cast<StringVertex*>(pVWEdge->getEnd());

			//std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
			//std::cout << pWVert->getID() << " w_edges: \n";
			EdgeDir transDir = !pVWEdge->getTwinDir();
			if(pWVert->getColor() == GC_GRAY)
			{
				EdgePtrVec w_edges = pWVert->getEdges(transDir);
				for(size_t j = 0; j < w_edges.size(); ++j)
				{
					//std::cout << "	edge: " << *w_edges[j] << "\n";
					StringEdge* pWXEdge = static_cast<StringEdge*>(w_edges[j]);
					size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
					if(trans_len <= longestLen)
					{
						if(pWXEdge->getEnd()->getColor() == GC_GRAY)
						{
							// X is the endpoint of an edge of V, therefore it is transitive
							pWXEdge->getEnd()->setColor(GC_BLACK);
							++marked_verts;
							//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive\n";
						}
					}
					else
						break;
				}
			}
		}

		// Stage 2
		for(size_t i = 0; i < edges.size(); ++i)
		{
			StringEdge* pVWEdge = static_cast<StringEdge*>(edges[i]);
			StringVertex* pWVert = static_cast<StringVertex*>(pVWEdge->getEnd());

			//std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
			//std::cout << pWVert->getID() << " w_edges: \n";
			EdgeDir transDir = !pVWEdge->getTwinDir();
			EdgePtrVec w_edges = pWVert->getEdges(transDir);
			for(size_t j = 0; j < w_edges.size(); ++j)
			{
				//std::cout << "	edge: " << *w_edges[j] << "\n";
				StringEdge* pWXEdge = static_cast<StringEdge*>(w_edges[j]);
				size_t len = pWXEdge->getSeqLen();

				if(len < FUZZ || j == 0)
				{
					if(pWXEdge->getEnd()->getColor() == GC_GRAY)
					{
						// X is the endpoint of an edge of V, therefore it is transitive
						pWXEdge->getEnd()->setColor(GC_BLACK);
						++marked_verts;
						//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive in stage 2\n";
					}
				}
				else
					break;
			}
		}

		bool trans_found = false;
		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_BLACK)
			{
				// Mark the edge and its twin for removal
				edges[i]->setColor(GC_BLACK);
				edges[i]->getTwin()->setColor(GC_BLACK);
				marked_edges += 2;
				trans_found = true;
			}
			edges[i]->getEnd()->setColor(GC_WHITE);
		}
	}

	return false;
}

// Remove all the marked edges
void SGTransRedVisitor::postvisit(StringGraph* pGraph)
{
	printf("TR marked %d verts and %d edges\n", marked_verts, marked_edges);
	pGraph->sweepEdges(GC_BLACK);
	assert(pGraph->checkColors(GC_WHITE));
}

void SGTrimVisitor::previsit(StringGraph* pGraph)
{
	num_island = 0;
	num_terminal = 0;
	num_contig = 0;
	pGraph->setColors(GC_WHITE);
}

// Mark any nodes that either dont have edges or edges in only one direction for removal
bool SGTrimVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	bool noext[2] = {0,0};

	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		if(pVertex->countEdges(dir) == 0)
		{
			pVertex->setColor(GC_BLACK);
			noext[idx] = 1;
		}
	}

	if(noext[0] && noext[1])
		num_island++;
	else if(noext[0] || noext[1])
		num_terminal++;
	else
		num_contig++;
	return noext[0] || noext[1];
}

// Remove all the marked edges
void SGTrimVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
	printf("island: %d terminal: %d contig: %d\n", num_island, num_terminal, num_contig);
}

//
// Bubble visit
//
void SGBubbleVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	num_bubbles = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGBubbleVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	bool bubble_found = false;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);
		if(edges.size() > 1)
		{
			// Check the vertices
			for(size_t i = 0; i < edges.size(); ++i)
			{
				Edge* pVWEdge = edges[i];
				Vertex* pWVert = pVWEdge->getEnd();

				// Get the edges from w in the same direction
				EdgeDir transDir = !pVWEdge->getTwinDir();
				EdgePtrVec wEdges = pWVert->getEdges(transDir);

				if(pWVert->getColor() == GC_RED)
					return false;

				// If the bubble has collapsed, there should only be one edge
				if(wEdges.size() == 1)
				{
					Vertex* pBubbleEnd = wEdges.front()->getEnd();
					if(pBubbleEnd->getColor() == GC_RED)
						return false;
				}
			}

			// Mark the vertices
			for(size_t i = 0; i < edges.size(); ++i)
			{
				Edge* pVWEdge = edges[i];
				Vertex* pWVert = pVWEdge->getEnd();

				// Get the edges from w in the same direction
				EdgeDir transDir = !pVWEdge->getTwinDir();
				EdgePtrVec wEdges = pWVert->getEdges(transDir);

				// If the bubble has collapsed, there should only be one edge
				if(wEdges.size() == 1)
				{
					Vertex* pBubbleEnd = wEdges.front()->getEnd();
					if(pBubbleEnd->getColor() == GC_BLACK)
					{
						// The endpoint has been visited, set this vertex as needing removal
						// and set the endpoint as unvisited
						pWVert->setColor(GC_RED);
						bubble_found = true;
					}
					else
					{
						// Endpoint has not been hit, set it to visited
						if(pBubbleEnd->getColor() == GC_RED || pWVert->getColor() == GC_RED)
						{
							std::cout << "Color stomp!\n";
						}
						pBubbleEnd->setColor(GC_BLACK);
						pWVert->setColor(GC_BLUE);
					}
				}
			}
			
			// Unmark vertices
			for(size_t i = 0; i < edges.size(); ++i)
			{
				Edge* pVWEdge = edges[i];
				Vertex* pWVert = pVWEdge->getEnd();

				// Get the edges from w in the same direction
				EdgeDir transDir = !pVWEdge->getTwinDir();
				EdgePtrVec wEdges = pWVert->getEdges(transDir);

				// If the bubble has collapsed, there should only be one edge
				if(wEdges.size() == 1)
				{
					Vertex* pBubbleEnd = wEdges.front()->getEnd();
					pBubbleEnd->setColor(GC_WHITE);
				}
				if(pWVert->getColor() == GC_BLUE)
					pWVert->setColor(GC_WHITE);
			}

			if(bubble_found)
				++num_bubbles;

		}
	}
	return bubble_found;
}

// Remove all the marked edges
void SGBubbleVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_RED);
	printf("bubbles: %d\n", num_bubbles);
	assert(pGraph->checkColors(GC_WHITE));
}

// Precondition: the graph has been transitively reduced
void SGVariantVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGVariantVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);

		// Perform a dependent pairwise comparison between nodes that do not share an edge
		for(size_t i = 0; i < edges.size(); ++i)
		{
			for(size_t j = i + 1; j < edges.size(); ++j)
			{
				const StringEdge* p_edgeI = CSE_CAST(edges[i]);
				const StringEdge* p_edgeJ = CSE_CAST(edges[j]);
				const StringVertex* p_vertI = CSV_CAST(p_edgeI->getEnd());
				const StringVertex* p_vertJ = CSV_CAST(p_edgeJ->getEnd());

				// Infer the edge comp
				EdgeComp comp = p_edgeI->getComp() == p_edgeJ->getComp() ? EC_SAME : EC_REVERSE;
				
				// Infer the i->j and j->i direction
				EdgeDir ij_dir = !p_edgeI->getTwin()->getDir();
				EdgeDir ji_dir = !p_edgeJ->getTwin()->getDir();

				// Ensure that there is not an edge between these already
				// hasEdge is not a particularly fast operation
				EdgeDesc ij_desc(p_edgeJ->getEndID(), ij_dir, comp);
				EdgeDesc ji_desc(p_edgeI->getEndID(), ji_dir, comp);
				if(p_edgeI->getEnd()->hasEdge(ij_desc) || p_edgeJ->getEnd()->hasEdge(ji_desc))
					continue;

				// Set up the matches between the root vertex and i/j
				// All coordinates are calculated from the point of view of pVertex
				Match match_i = p_edgeI->getMatch();
				Match match_j = p_edgeJ->getMatch();
	
				// Infer the match_ij based match_i and match_j
				Match match_ij = Match::infer(match_i, match_j);
	
				std::string seq_i = p_vertI->getSeq();
				std::string seq_j = !match_ij.isRC() ? p_vertJ->getSeq() : reverseComplement(p_vertJ->getSeq());
				
				// Expand the match outwards so that one sequence is left terminal 
				// and one is right terminal
				match_ij.canonize();
				match_ij.expand();
				
				// Extract the match strings
				std::string ms_i = match_ij.coord[0].getSubstring(seq_i);
				std::string ms_j = match_ij.coord[1].getSubstring(seq_j);

				std::cout << "COMP: " << comp << "\n";
				std::cout << "match_i: " << match_i << " " << p_edgeI->getComp() << "\n";
				std::cout << "match_j: " << match_j << " " << p_edgeJ->getComp() << "\n";
				std::cout << "match_ij: " << match_ij << "\n";
				std::cout << "ms_i: " << ms_i << "\n";
				std::cout << "ms_j: " << ms_j << "\n";
				/*
				Overlap o = inferOverlap(seq_i, seq_j, start_i, start_j, p_vertI->getID(), p_vertJ->getID(), comp == EC_REVERSE);
				m_fileHandle << o << "\n";
				*/
			}
		}
	}
	return false;
}

// Simple gapless count of mismatches between seq_i and seq_j given the starting position of the matches
Overlap SGVariantVisitor::inferOverlap(const std::string& seq_i, const std::string& seq_j, int start_i, int start_j, 
                                    const std::string& id_i, const std::string& id_j, bool isReverse) const
{
	int match_span = std::min(seq_i.length() - start_i, seq_j.length() - start_j);
	int num_diff = 0;
	printf("seq_i(%d): %s\nseq_j(%d): %s\n", start_i, seq_i.c_str(), start_j, seq_j.c_str());
	for(int k = 0; k < match_span; ++k)
	{
		if(seq_i[start_i + k] != seq_j[start_j + k])
			++num_diff;
		printf("%c - %c\n", seq_i[start_i + k], seq_j[start_j + k]);
	}

	// Create the overlap
	SeqCoord sc_i(seq_i.length(), start_i, start_i + match_span - 1);
	SeqCoord sc_j(seq_j.length(), start_j, start_j + match_span - 1);
	Overlap o(id_i, sc_i, id_j, sc_j, isReverse, num_diff);
	return o;
}

//
void SGVariantVisitor::postvisit(StringGraph* /*pGraph*/)
{
}

void SGErrorRemovalVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);

}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGErrorRemovalVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	EdgePtrVec edges = pVertex->getEdges();
	int bp_overlap = 0;
	int num_diff = 0;
	for(size_t i = 0; i < edges.size(); ++i)
	{
		const StringEdge* pSE = CSE_CAST(edges[i]);
		bp_overlap += pSE->getMatchLength();
		num_diff += pSE->getNumDiff();
	}

	double ner = (double)num_diff / bp_overlap;
	//printf("ND: %d NO: %d NER: %lf\n", num_diff, bp_overlap, ner);

	if(ner > m_maxErrorRate)
		pVertex->setColor(GC_BLACK);

	return false;
}

// Remove all the marked edges
void SGErrorRemovalVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
}


void SGGraphStatsVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_terminal = 0;
	num_island = 0;
	num_monobranch = 0;
	num_dibranch = 0;
	num_transitive = 0;
	num_edges = 0;
	num_vertex = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGGraphStatsVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	int s_count = pVertex->countEdges(ED_SENSE);
	int as_count = pVertex->countEdges(ED_ANTISENSE);
	if(s_count == 0 && as_count == 0)
		++num_island;
	else if(s_count == 0 || as_count == 0)
		++num_terminal;

	if(s_count > 1 && as_count > 1)
		++num_dibranch;
	else if(s_count > 1 || as_count > 1)
		++num_monobranch;

	if(s_count == 1 || as_count == 1)
		++num_transitive;

	num_edges += (s_count + as_count);
	++num_vertex;
	return false;
}

// Remove all the marked edges
void SGGraphStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("island: %d terminal: %d monobranch: %d dibranch: %d transitive: %d\n", num_island, num_terminal,
	                                                                               num_monobranch, num_dibranch, num_transitive);
	printf("Total Vertices: %d Total Edges: %d\n", num_vertex, num_edges);
}

