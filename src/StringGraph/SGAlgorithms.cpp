//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#include "SGAlgorithms.h"
#include "SGUtil.h"

// Visitor which outputs the graph in fasta format
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	m_fileHandle << ">" << pVertex->getID() << " " <<  pVertex->getSeq().length() 
				 << " " << pVertex->getReadCount() << "\n";
	m_fileHandle << pVertex->getSeq() << "\n";
	return false;
}

void SGTransRedVisitor::previsit(StringGraph* pGraph)
{
	// Set all the vertices in the graph to "vacant"
	pGraph->setColors(GC_WHITE);
	pGraph->sortVertexAdjListsByLen();

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
		EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted

		if(edges.size() == 0)
			continue;

		for(size_t i = 0; i < edges.size(); ++i)
			(edges[i])->getEnd()->setColor(GC_GRAY);

		Edge* pLongestEdge = edges.back();
		size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;
		
		// Stage 1
		for(size_t i = 0; i < edges.size(); ++i)
		{
			Edge* pVWEdge = edges[i];
			Vertex* pWVert = pVWEdge->getEnd();

			//std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
			//std::cout << pWVert->getID() << " w_edges: \n";
			EdgeDir transDir = !pVWEdge->getTwinDir();
			if(pWVert->getColor() == GC_GRAY)
			{
				EdgePtrVec w_edges = pWVert->getEdges(transDir);
				for(size_t j = 0; j < w_edges.size(); ++j)
				{
					Edge* pWXEdge = w_edges[j];
					size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
					if(trans_len <= longestLen)
					{
						if(pWXEdge->getEnd()->getColor() == GC_GRAY)
						{
							// X is the endpoint of an edge of V, therefore it is transitive
							pWXEdge->getEnd()->setColor(GC_BLACK);
							++marked_verts;
							//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << "\n";
						}
					}
					else
						break;
				}
			}
		}
		
		WARN_ONCE("Transitive Reduction stage 2 disabled")
		/*
		// Stage 2
		for(size_t i = 0; i < edges.size(); ++i)
		{
			Edge* pVWEdge = edges[i];
			Vertex* pWVert = pVWEdge->getEnd();

			//std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
			//std::cout << pWVert->getID() << " w_edges: \n";
			EdgeDir transDir = !pVWEdge->getTwinDir();
			EdgePtrVec w_edges = pWVert->getEdges(transDir);
			for(size_t j = 0; j < w_edges.size(); ++j)
			{
				//std::cout << "	edge: " << *w_edges[j] << "\n";
				Edge* pWXEdge = w_edges[j];
				size_t len = pWXEdge->getSeqLen();

				if(len < FUZZ || j == 0)
				{
					if(pWXEdge->getEnd()->getColor() == GC_GRAY)
					{
						// X is the endpoint of an edge of V, therefore it is transitive
						pWXEdge->getEnd()->setColor(GC_BLACK);
						++marked_verts;
						//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << " in stage 2\n";
					}
				}
				else
					break;
			}
		}
		*/

		bool trans_found = false;
		size_t trans_count = 0;
		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_BLACK)
			{
				// Mark the edge and its twin for removal
				edges[i]->setColor(GC_BLACK);
				edges[i]->getTwin()->setColor(GC_BLACK);
				marked_edges += 2;
				trans_found = true;
				trans_count++;
			}
			edges[i]->getEnd()->setColor(GC_WHITE);
		}
		/*
		if(trans_count + 1 != edges.size())
		{
			printf("Vertex %s could not be completely reduced (%d, %d)\n", pVertex->getID().c_str(), (int)trans_count, (int)edges.size());
			for(size_t i = 0; i < edges.size(); ++i)
			{
				if(edges[i]->getColor() != GC_BLACK)
					std::cout << "Remaining edge: " << *edges[i] << "\n";
			}
		}
		*/
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
			//std::cout << "Found terminal: " << pVertex->getID() << "\n";
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

// Detect and remove duplicate edges
bool SGDuplicateVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	pVertex->makeUnique();
	return false;
}

void SGIslandVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

// Mark any nodes that dont have edges
bool SGIslandVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	if(pVertex->countEdges() == 0)
	{
		pVertex->setColor(GC_BLACK);
		return true;
	}
	return false;
}

// Remove all the marked vertices
void SGIslandVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
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

// 
void SGVariantVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

// 
bool SGVariantVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	// Skip vertex if its been marked as contained
	if(pVertex->getColor() == GC_RED)
		return false;
	bool changed_graph = false;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);

		// Perform a dependent pairwise comparison between nodes that do not share an edge
		for(size_t i = 0; i < edges.size(); ++i)
		{
			for(size_t j = i + 1; j < edges.size(); ++j)
			{
				const Edge* p_edgeI = edges[i];
				const Edge* p_edgeJ = edges[j];

				/*
				if(p_edgeI->getColor() == GC_BLACK || p_edgeJ->getColor() == GC_BLACK)
				{
					//std::cout << "Skipping newly created edges\n";
					continue;
				}
				*/
				const Vertex* p_vertI = p_edgeI->getEnd();
				const Vertex* p_vertJ = p_edgeJ->getEnd();

				// Infer the edge comp
				EdgeComp comp = (p_edgeI->getComp() == p_edgeJ->getComp()) ? EC_SAME : EC_REVERSE;
				
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
				std::string seq_j = p_vertJ->getSeq();
				
				// Expand the match outwards so that one sequence is left terminal 
				// and one is right terminal
				match_ij.expand();
				
				// Extract the match strings
				std::string ms_i = match_ij.coord[0].getSubstring(seq_i);
				std::string ms_j = match_ij.coord[1].getSubstring(seq_j);

				if(match_ij.isRC())
					ms_j = reverseComplement(ms_j);
				
				/*
				std::cout << "COMP: " << comp << "\n";
				std::cout << "match_i: " << match_i << " " << p_edgeI->getComp() << "\n";
				std::cout << "match_j: " << match_j << " " << p_edgeJ->getComp() << "\n";
				std::cout << "match_ij: " << match_ij << "\n";
				std::cout << "ms_i: " << ms_i << "\n";
				std::cout << "ms_j: " << ms_j << "\n";
				*/
				int numDiffs = countDifferences(ms_i, ms_j, ms_i.length());
				match_ij.setNumDiffs(numDiffs);

				//std::cout << "diffs: " << match_ij.getNumDiffs() << "\n";
				double error_rate = (double)numDiffs / ms_i.length();
				if(error_rate < 0.1)
				{
					const std::string& id_i = p_edgeI->getEndID();
					const std::string& id_j = p_edgeJ->getEndID();
					Overlap o(id_i, id_j, match_ij);

					// Ensure there isnt a containment relationship
					if(o.match.coord[0].isContained() || o.match.coord[1].isContained())
					{
						// Mark the contained vertex for deletion
						size_t idx = o.getContainedIdx();

						/*
						std::cout << "ci: " << seq_i << "\n";
						std::cout << "cj: " << seq_j << "\n";
						std::cout << "ms_i: " << ms_i << "\n";
						std::cout << "ms_j: " << ms_j << "\n";
						*/
						if(idx == 0)
							p_edgeI->getEnd()->setColor(GC_RED);
						else
							p_edgeJ->getEnd()->setColor(GC_RED);
						
					}
					else
					{
						// add the edges to the graph
						Edge* p_edgeIJ = createEdges(pGraph, o);
						Edge* p_edgeJI = p_edgeIJ->getTwin();
						p_edgeIJ->setColor(GC_BLACK);
						p_edgeJI->setColor(GC_BLACK);
						//std::cout << "Created edge " << *p_edgeIJ << "\n";
						//std::cout << "Created edge " << *p_edgeJI << "\n";
						assert(p_edgeIJ->getComp() == comp);
					}
					changed_graph = true;
				}
			}
		}
	}
	return changed_graph;
}

//
void SGVariantVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_RED);
}

// 
void SGTCVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	SGEdgeClassVisitor ecv;
	pGraph->visit(ecv);
	
	ngb = ecv.getNumGood();
	nbb = ecv.getNumBad();
}

// 
bool SGTCVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	// Skip vertex if its been marked as contained
	if(pVertex->getColor() == GC_RED)
		return false;

	bool changed_graph = false;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);

		// Perform a dependent pairwise comparison between nodes that do not share an edge
		for(size_t i = 0; i < edges.size(); ++i)
		{
			for(size_t j = i + 1; j < edges.size(); ++j)
			{
				const Edge* p_edgeI = edges[i];
				const Edge* p_edgeJ = edges[j];
				const Vertex* p_vertI = p_edgeI->getEnd();
				const Vertex* p_vertJ = p_edgeJ->getEnd();

				// Infer the edge comp
				EdgeComp comp = (p_edgeI->getComp() == p_edgeJ->getComp()) ? EC_SAME : EC_REVERSE;
				
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
				std::string seq_j = p_vertJ->getSeq();
				
				// Expand the match outwards so that one sequence is left terminal 
				// and one is right terminal
				match_ij.expand();
				
				// Extract the unmatched portion of each strings
				std::string cs_i = match_i.coord[1].getComplementString(seq_i);
				std::string cs_j = match_j.coord[1].getComplementString(seq_j);

				if(match_ij.isRC())
					cs_j = reverseComplement(cs_j);
				
				if(match_i.coord[0].isLeftExtreme())
				{
					cs_i = reverse(cs_i);
					cs_j = reverse(cs_j);
				}

				int min_size = std::min(cs_i.length(), cs_j.length());
				
				/*
				std::cout << "COMP: " << comp << "\n";
				std::cout << "match_i: " << match_i << " " << p_edgeI->getComp() << "\n";
				std::cout << "match_j: " << match_j << " " << p_edgeJ->getComp() << "\n";
				std::cout << "match_ij: " << match_ij << "\n";
				std::cout << "cs_i: " << cs_i << "\n";
				std::cout << "cs_j: " << cs_j << "\n";
				(void)pGraph;
				*/
				int numDiffs = countDifferences(cs_i, cs_j, min_size);
				match_ij.setNumDiffs(numDiffs);

				//std::cout << "diffs: " << match_ij.getNumDiffs() << "\n";
				//double error_rate = (double)numDiffs / min_size;
				if(numDiffs <= 2)
				{
					const std::string& id_i = p_edgeI->getEndID();
					const std::string& id_j = p_edgeJ->getEndID();
					Overlap o(id_i, id_j, match_ij);

					// Ensure there isnt a containment relationship
					if(o.match.coord[0].isContained() || o.match.coord[1].isContained())
					{
						// Mark the contained vertex for deletion
						size_t idx = o.getContainedIdx();
						/*
						std::cout << "ci: " << seq_i << "\n";
						std::cout << "cj: " << seq_j << "\n";
						std::cout << "ms_i: " << ms_i << "\n";
						std::cout << "ms_j: " << ms_j << "\n";
						*/
						if(idx == 0)
							p_edgeI->getEnd()->setColor(GC_RED);
						else
							p_edgeJ->getEnd()->setColor(GC_RED);
						
					}
					else
					{
						// add the edges to the graph
						Edge* p_edgeIJ = createEdges(pGraph, o);
						Edge* p_edgeJI = p_edgeIJ->getTwin();
						p_edgeIJ->setColor(GC_BLACK);
						p_edgeJI->setColor(GC_BLACK);
						//std::cout << "Created edge " << *p_edgeIJ << "\n";
						//std::cout << "Created edge " << *p_edgeJI << "\n";
						assert(p_edgeIJ->getComp() == comp);
					}
					changed_graph = true;
				}
			}
		}
	}
	return changed_graph;
}

//
void SGTCVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_RED);

	SGEdgeClassVisitor ecv;
	pGraph->visit(ecv);
	
	int nga = ecv.getNumGood();
	int nba = ecv.getNumBad();

	int ngc = nga - ngb;
	int nbc = nba - nbb;
	double rel = (double)ngc / (nbc + ngc);

	printf("ng before: %d ng after: %d (%d) nb before: %d nb after: %d (%d) prop good: %lf\n", ngb, nga, ngc, nbb, nba, nbc, rel);
}

void SGVertexPairingVisitor::previsit(StringGraph* pGraph)
{
	num_paired = 0;
	num_unpaired = 0;
	pGraph->setColors(GC_WHITE);
}

// Visit each vertex in the graph, find its pair and link them
bool SGVertexPairingVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	// do nothing if the pairing was already set
	if(pVertex->getPairVertex() == NULL)
	{
		std::string id = pVertex->getID();
		std::string pid = getPairID(id);
		Vertex* pPairV = pGraph->getVertex(pid);
		if(pPairV != NULL)
		{
			pVertex->setPairVertex(pPairV);
			pPairV->setPairVertex(pVertex);
			num_paired++;
		}
		else
		{
			//pVertex->setColor(GC_BLACK);
			num_unpaired++;
		}
	}
	return false;
}

void SGVertexPairingVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("Graph has %d paired vertices, %d verts with no pairs\n", num_paired, num_unpaired);
	//pGraph->sweepVertices(GC_BLACK);
}

//
void SGPETrustVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	//printf("TOKEN\tTRUSTED\tNOT\tDIFFSTRAND\tTOTAL\n");
}


// Visit each vertex in the graph and determine which edges are supported through
// read pairing
bool SGPETrustVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	Vertex* pPairVertex = pVertex->getPairVertex();
	if(pPairVertex == NULL)
		return false;

	// First, mark all pair vertices that overlap the pair of this node
	// The set of marked vertices that overlap pVertex are the trusted vertices
	EdgePtrVec pairEdgeVec = pPairVertex->getEdges();
	for(size_t i = 0; i < pairEdgeVec.size(); ++i)
	{
		// Get the pair of the endpoint of this edge
		Vertex* pBackVertex = pairEdgeVec[i]->getEnd()->getPairVertex();
		if(pBackVertex != NULL)
			pBackVertex->setColor(GC_RED);
	}

	EdgePtrVec vertEdgeVec = pVertex->getEdges();
	
	bool changed = true;
	while(changed)
	{
		changed = false;
		// Propogate trust
		for(size_t i = 0; i < vertEdgeVec.size(); ++i)
		{
			Vertex* pCurr = vertEdgeVec[i]->getEnd();
			if(pCurr->getColor() != GC_RED)
			{
				// If any vertex that pCurr overlaps with is red, mark it red too
				EdgePtrVec currEdgeVec = pCurr->getEdges();
				for(size_t j = 0; j < currEdgeVec.size(); ++j)
				{
					if(currEdgeVec[j]->getEnd()->getColor() == GC_RED)
					{
						pCurr->setColor(GC_RED);
						changed = true;
						break;
					}
				}
			}
		}
	}

	// 
	int trusted = 0;
	int nottrusted = 0;
	int diffstrand = 0;
	for(size_t i = 0; i < vertEdgeVec.size(); ++i)
	{
		if(vertEdgeVec[i]->getEnd()->getColor() == GC_RED)
		{
			trusted++;
			vertEdgeVec[i]->isTrusted = true;
		}
		else
		{
			nottrusted++;
		}
	}

	(void)diffstrand;
	//printf("TOKEN\t%d\t%d\t%d\t%zu\n", trusted, nottrusted, diffstrand, vertEdgeVec.size());

	// Reset all the vertex colors
	for(size_t i = 0; i < pairEdgeVec.size(); ++i)
	{
		// Get the pair of the endpoint of this edge
		Vertex* pBackVertex = pairEdgeVec[i]->getEnd()->getPairVertex();
		if(pBackVertex)
			pBackVertex->setColor(GC_WHITE);
	}

	for(size_t i = 0; i < vertEdgeVec.size(); ++i)
		vertEdgeVec[i]->getEnd()->setColor(GC_WHITE);
	
	return false;
}

//
void SGPETrustVisitor::postvisit(StringGraph* /*pGraph*/)
{

}

bool SGPairedOverlapVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	Vertex* pPairSV = pVertex->getPairVertex();
	if(pPairSV == NULL)
		return false;

	EdgePtrVec edges = pVertex->getEdges();

	// Determine which vertices that are paired to pVertex
	// have a pair that overlaps with pPairVertex
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pVWEdge = edges[i];
		Vertex* pW = pVWEdge->getEnd();
		Vertex* pPairW = pW->getPairVertex();
		if(pPairW == NULL)
			continue;

		EdgePtrVec ppw_edges = pPairW->findEdgesTo(pPairSV->getID());
		size_t overlap_len = pVWEdge->getMatchLength();

		if(pVWEdge->getComp() == EC_SAME)
		{
			if(ppw_edges.size() == 1)
			{
				Edge* pPPEdge = ppw_edges.front();
				size_t pair_overlap_len = pPPEdge->getMatchLength();
				printf("pairoverlap\t%s\t%s\t%zu\t%zu\n", pVertex->getID().c_str(), pW->getID().c_str(), overlap_len, pair_overlap_len);
			}
			else
			{
				printf("pairoverlap\t%s\t%s\t%zu\t%d\n", pVertex->getID().c_str(), pW->getID().c_str(), overlap_len, 0);
			}
		}
	}
	return false;
}

void SGPEResolveVisitor::previsit(StringGraph*)
{
	printf("RESOLVE\tDIR\tPOS\tGOOD_RED\tGOOD_BLACK\tBAD_RED\tBAD_BLACK\tTOTAL\tCONFLICTED\tRESOLVABLE\n");
}

bool SGPEResolveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);

		int num_bad_black = 0;
		int num_good_black = 0;
		int num_bad_red = 0;
		int num_good_red = 0;
		
		for(size_t i = 0; i < edges.size(); ++i)
		{
			Edge* pVWEdge = edges[i];
			Vertex* pW = pVWEdge->getEnd();

			// Get the distance between the vertices using the debug positions
			int distance = abs(pW->dbg_position - pVertex->dbg_position);
			int overlap_len = pVWEdge->getMatchLength();
			int sum = distance + overlap_len;
			
			bool isGood = false;
			if(sum == (int)pVertex->getSeq().size())
			{
				if(!pVWEdge->isTrusted)
					num_good_black++;
				else
					num_good_red++;
				isGood = true;
			}
			else
			{
				if(pVWEdge->isTrusted)
					num_bad_black++;
				else
					num_bad_red++;			
				isGood = false;
			}
		}

		bool isConflicted = edges.size() > 1;
		bool isResolvable = isConflicted && num_good_red > 0 && num_bad_red == 0;

		printf("RESOLVE\t%zu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", idx, pVertex->dbg_position, num_good_red, num_good_black, num_bad_red, num_bad_black, (int)edges.size(), isConflicted, isResolvable);
	}
	return false;
}

void SGPEConflictRemover::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	num_same = 0;
	num_diff = 0;
}

bool SGPEConflictRemover::visit(StringGraph* pGraph, Vertex* pVertex)
{
	(void)pGraph;
	(void)pVertex;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		EdgePtrVec edges = pVertex->getEdges(dir);

		if(edges.size() > 1)
		{
			bool hasTrusted = false;
			for(size_t j = 0; j < edges.size(); ++j)
			{
				if(edges[j]->isTrusted)
				{
					hasTrusted = true;
				}
			}
			
			if(hasTrusted)
			{
				for(size_t j = 0; j < edges.size(); ++j)
				{
					if(!edges[j]->isTrusted)
					{
						edges[j]->setColor(GC_BLACK);
						edges[j]->getTwin()->setColor(GC_BLACK);
					}
					if(edges[j]->getComp() == EC_SAME)
						num_same++;
					else
						num_diff++;
				}
			}
		}	
	}
	return 0;
}

void SGPEConflictRemover::postvisit(StringGraph* pGraph)
{
	pGraph->sweepEdges(GC_BLACK);
	printf("Removed %d diff %d same\n", num_diff, num_same);
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


void SGEdgeClassVisitor::previsit(StringGraph* /*pGraph*/)
{
	num_good = 0;
	num_bad = 0;
	num_conflicted = 0;
	num_trusted = 0;
	num_nottrusted = 0;
}

bool SGEdgeClassVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
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
void SGEdgeClassVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("Num good: %d Num bad: %d Num conflicted: %d Num trusted: %d Num not trusted: %d\n", num_good, num_bad, num_conflicted, num_trusted, num_nottrusted);
}
