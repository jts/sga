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

// Infer an overlap from two edges
// The input edges are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap SGAlgorithms::inferTransitiveOverlap(const Edge* pXY, const Edge* pYZ)
{
	// Construct the match
	Match match_yx = pXY->getMatch();
	match_yx.swap(); 
	Match match_yz = pYZ->getMatch();

	// Infer the match_ij based match_i and match_j
	Match match_xz = Match::infer(match_yx, match_yz);
	match_xz.expand();

	// Convert the match to an overlap
	Overlap ovr(pXY->getStartID(), pYZ->getEndID(), match_xz);
	return ovr;
}

// Return true if XZ has an overlap
bool SGAlgorithms::hasTransitiveOverlap(const Edge* pXY, const Edge* pYZ)
{
	// ensure the y component is the first one in the seqcoord
	// as the match algorithms expect the common basis to be the first
	// component
	Match match_yx = pXY->getMatch();
	match_yx.swap(); 
	Match match_yz = pYZ->getMatch();
	return Match::doMatchesIntersect(match_yx, match_yz);
}

//
// SGFastaVisitor - output the vertices in the graph in 
// fasta format
//
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	m_fileHandle << ">" << pVertex->getID() << " " <<  pVertex->getSeq().length() 
				 << " " << pVertex->getReadCount() << "\n";
	m_fileHandle << pVertex->getSeq() << "\n";
	return false;
}


//
// SGOverlapWriterVisitor - write all the overlaps in the graph to a file 
//
bool SGOverlapWriterVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Overlap ovr = edges[i]->getOverlap();
		if(ovr.id[0] < ovr.id[1])
			m_fileHandle << ovr << "\n";
	}
	return false;
}




//
// SGTransRedVisitor - Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
void SGTransRedVisitor::previsit(StringGraph* pGraph)
{
	// The graph must not have containments
	assert(!pGraph->hasContainment());

	// Set all the vertices in the graph to "vacant"
	pGraph->setColors(GC_WHITE);
	pGraph->sortVertexAdjListsByLen();

	marked_verts = 0;
	marked_edges = 0;
}

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

//
// SGContainRemoveVisitor - Removes contained
// vertices from the graph
//
void SGContainRemoveVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

//
bool SGContainRemoveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Match m = edges[i]->getMatch();
		if(m.isContainment())
		{
			if(pVertex->getID() > edges[i]->getEnd()->getID())
			{
				pVertex->setColor(GC_BLACK);
				break;
			}
		}
	}
	return false;
}

void SGContainRemoveVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
	pGraph->setContainmentFlag(false);
}
//
// SGRealignVisitor - Infer locations of potentially
// missing edges and add them to the graph
//
void SGRealignVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

bool SGRealignVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	bool graph_changed = false;
	const int MIN_OVERLAP = 41;
	const double MAX_ERROR = 0.20;
	static int visited = 0;
	++visited;
	if(visited % 50000 == 0)
		std::cout << "visited: " << visited << "\n";

	// Explore the neighborhood around this graph for potentially missing overlaps
	CandidateVector candidates = getMissingCandidates(pGraph, pVertex, MIN_OVERLAP);

	for(size_t i = 0; i < candidates.size(); ++i)
	{
		Candidate& c = candidates[i];
		int numDiff = c.ovr.match.countDifferences(pVertex->getSeq(), c.pEndpoint->getSeq());
		double error_rate = double(numDiff) / double(c.ovr.match.getMinOverlapLength());

		//std::cout << "Actual:\n";
		//c.ovr.match.printMatch(pVertex->getSeq(), c.pEndpoint->getSeq());
		//std::cout << "ER: " << error_rate << "\n";
		//std::cout << "OVR:\t" << c.ovr << "\n";
		
		if(error_rate < MAX_ERROR)
		{
			Edge* p_edgeXZ = createEdges(pGraph, c.ovr, true);
			Edge* p_edgeZX = p_edgeXZ->getTwin();
			p_edgeXZ->setColor(GC_WHITE);
			p_edgeZX->setColor(GC_WHITE);
			graph_changed = true;
		}
		/*else
		{
			std::cout << "Rejected actual:\n";
			c.ovr.match.printMatch(pVertex->getSeq(), c.pEndpoint->getSeq());
			std::cout << "ER: " << error_rate << "\n";
			std::cout << "OVR:\t" << c.ovr << "\n";
		}*/
	}
	return graph_changed;
}

// Explore the neighborhood around a vertex looking for missing overlaps
SGRealignVisitor::CandidateVector SGRealignVisitor::getMissingCandidates(StringGraph* /*pGraph*/, Vertex* pVertex, int minOverlap) const
{
	CandidateVector out;

	// Mark the vertices that are reached from this vertex as black to indicate
	// they already are overlapping
	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	{
		edges[i]->getEnd()->setColor(GC_BLACK);
	}
	pVertex->setColor(GC_BLACK);

	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pXY = edges[i];
		EdgePtrVec neighborEdges = pXY->getEnd()->getEdges();
		for(size_t j = 0; j < neighborEdges.size(); ++j)
		{
			Edge* pYZ = neighborEdges[j];
			if(pYZ->getEnd()->getColor() != GC_BLACK)
			{
				// Infer the overlap object from the edges
				if(SGAlgorithms::hasTransitiveOverlap(pXY, pYZ))
				{
					Overlap ovr_xz = SGAlgorithms::inferTransitiveOverlap(pXY, pYZ);
					if(ovr_xz.match.getMinOverlapLength() >= minOverlap)
					{
						out.push_back(Candidate(pYZ->getEnd(), ovr_xz));
						pYZ->getEnd()->setColor(GC_BLACK);
					}
				}
			}
		}
	}

	// Reset colors
	for(size_t i = 0; i < edges.size(); ++i)
		edges[i]->getEnd()->setColor(GC_WHITE);
	pVertex->setColor(GC_WHITE);
	for(size_t i = 0; i < out.size(); ++i)
		out[i].pEndpoint->setColor(GC_WHITE);
	return out;
}

//
void SGRealignVisitor::postvisit(StringGraph*)
{
}


//
// SGCloseGroupVisitor - Attempt to close transitive groups
// by inferring missing edges from the graph that result from a 
// high error rate.
//
void SGGroupCloseVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	numGroupsOpen = 0;
	numGroupsClosed = 0;
	numEdgesRejected = 0;
}

bool SGGroupCloseVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	(void)pGraph;
	bool changed = false;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		TransitiveGroupCollection tgc = pVertex->computeTransitiveGroups(dir);
		if(tgc.numGroups() <= 1)
		{
			++numGroupsClosed;
			//continue;
		}
		else
		{
			++numGroupsOpen;
		}

		EdgePtrVec edges = pVertex->getEdges(dir);
/*
		TransitiveGroup& headGroup = tgc.getGroup(0);
		for(size_t j = 1; j < tgc.numGroups(); ++j)
		{
			TransitiveGroup& secondGroup = tgc.getGroup(j);

			WARN_ONCE("SGGroupCloseVisitor -- DIFF GROUPS and HANDLE CONTAINMENTS");

			Edge* pCandidate = secondGroup.getIrreducible();
			std::string candCons = pCandidate->getEnd()->getInferredConsensus();

			for(size_t i = 0; i < headGroup.numElements(); ++i)
			{
				Edge* pTarget = headGroup.getEdge(i);
*/
		// Perform a dependent pairwise comparison between nodes that do not share an edge
		for(size_t i = 0; i < edges.size(); ++i)
		{
			for(size_t j = i + 1; j < edges.size(); ++j)
			{
				Edge* pCandidate = edges[i];
				Edge* pTarget = edges[j];

				// Ensure these vertices dont already have an edge
				// Infer the edge comp
				EdgeComp comp = (pCandidate->getComp() == pTarget->getComp()) ? EC_SAME : EC_REVERSE;
				
				// Infer the i->j and j->i direction
				EdgeDir ij_dir = !pCandidate->getTwin()->getDir();
				EdgeDir ji_dir = !pTarget->getTwin()->getDir();

				// Ensure that there is not an edge between these already
				// hasEdge is not a particularly fast operation
				EdgeDesc ij_desc(pTarget->getEndID(), ij_dir, comp);
				EdgeDesc ji_desc(pCandidate->getEndID(), ji_dir, comp);

				if(pCandidate->getEnd()->hasEdge(ij_desc) || pTarget->getEnd()->hasEdge(ji_desc))
					continue;

				std::string candCons = pCandidate->getEnd()->getInferredConsensus();
				std::string targetCons = pTarget->getEnd()->getInferredConsensus();
				
				// Construct the match
				Match match_i = pCandidate->getMatch();
				Match match_j = pTarget->getMatch();

				// Infer the match_ij based match_i and match_j
				Match match_ij = Match::infer(match_i, match_j);
				match_ij.expand();

				// Convert the match to an overlap
				Overlap ovr(pCandidate->getEndID(), pTarget->getEndID(), match_ij);
				//std::cout << "CANDCONS: " << candCons << "\n";
				//std::cout << "TARGCONS: " << targetCons << "\n";

				//
				int numDiffCons = ovr.match.countDifferences(candCons, targetCons);
				bool accepted = (numDiffCons < 3);
				/*
				if(numDiffCons > 5)
				{
					std::cout << "Consensus:\n";
					ovr.match.printMatch(candCons, targetCons);
				
					std::cout << "Actual:\n";
					ovr.match.printMatch(pCandidate->getEnd()->getSeq(), pTarget->getEnd()->getSeq());
					
					int numDiffActual = ovr.match.countDifferences(pCandidate->getEnd()->getSeq(), pTarget->getEnd()->getSeq());
					printf("GCV\t%d\t%d\t%d\n", numDiffActual, numDiffCons, accepted);
				}
				*/
				//printf("GCV\t%d\t%d\t%d\n", numDiffActual, numDiffCons, accepted);
				
				if(accepted)
				{
					if(ovr.match.isContainment())
					{
						// Mark edge and vertex for removal
						//pCandidate->setColor(GC_BLACK);
						//pCandidate->getEnd()->setColor(GC_BLACK);
					}
					else
					{
						Edge* p_edgeIJ = createEdges(pGraph, ovr);
						Edge* p_edgeJI = p_edgeIJ->getTwin();
						p_edgeIJ->setColor(GC_WHITE);
						p_edgeJI->setColor(GC_WHITE);
						changed = true;
					}
				}
				else
				{
					++numEdgesRejected;
				}
			}
		}
	}
	return changed;
}

// Remove all the marked edges
void SGGroupCloseVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepVertices(GC_BLACK);
	printf("Num groups closed: %d Num groups open: %d Num edges rejected: %d\n", numGroupsClosed, numGroupsOpen, numEdgesRejected);
}

//
// SGTrimVisitor - Remove "dead-end" vertices from the graph
//
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

//
// SGDuplicateVisitor - Detect and remove duplicate edges
//
bool SGDuplicateVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	pVertex->makeUnique();
	return false;
}

//
// SGIslandVisitor - Remove island (unconnected) vertices
//
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
// SGBubbleVisitor - Find and collapse variant
// "bubbles" in the graph
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
// SGTCVisitor - Perform a transitive closure on the graph, 
// inferring edges that are missing because of inexact overlaps
// 
void SGTCVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
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
	
				// Expand the match outwards so that one sequence is left terminal 
				// and one is right terminal
				match_ij.expand();

				std::string seq_i = p_vertI->getSeq();
				std::string seq_j = p_vertJ->getSeq();
				

				
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
}


//
// SGGraphStatsVisitor - Collect summary stasitics
// about the graph
//
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

