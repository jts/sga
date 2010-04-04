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
#include "ErrorCorrect.h"

// Infer an overlap from two edges
// The input edges are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap SGAlgorithms::inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
	// Construct the match
	Match match_yx = ovrXY.match;
	match_yx.swap(); 
	Match match_yz = ovrYZ.match;

	// Infer the match_ij based match_i and match_j
	Match match_xz = Match::infer(match_yx, match_yz);
	match_xz.expand();

	// Convert the match to an overlap
	Overlap ovr(ovrXY.id[0], ovrYZ.id[1], match_xz);
	return ovr;
}

// Return true if XZ has an overlap
bool SGAlgorithms::hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
	// ensure the y component is the first one in the seqcoord
	// as the match algorithms expect the common basis to be the first
	// component
	Match match_yx = ovrXY.match;
	match_yx.swap(); 
	Match match_yz = ovrYZ.match;
	return Match::doMatchesIntersect(match_yx, match_yz);
}

// Construct an extended multioverlap for a vertex
MultiOverlap SGAlgorithms::makeExtendedMultiOverlap(const Vertex* pVertex)
{
	VertexOverlapSet overlapSet;
	findOverlapSet(pVertex, overlapSet);

	MultiOverlap mo(pVertex->getID(), pVertex->getSeq());
	for(VertexOverlapSet::const_iterator iter = overlapSet.begin();
	    iter != overlapSet.end(); ++iter)
	{
		mo.add(iter->pVertex->getSeq(), iter->ovr);
	}
	return mo;
}

//
void SGAlgorithms::makeExtendedSeqTries(const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie)
{
	double lp = log(p_error);
	VertexOverlapSet overlapSet;
	findOverlapSet(pVertex, overlapSet);

	for(VertexOverlapSet::const_iterator iter = overlapSet.begin();
	    iter != overlapSet.end(); ++iter)
	{
		// Coord[0] of the match is wrt pVertex, coord[1] is the other read
		std::string overlapped = iter->ovr.match.coord[1].getSubstring(iter->pVertex->getSeq());
		if(iter->ovr.match.isRC())
			overlapped = reverseComplement(overlapped);

		if(iter->ovr.match.coord[0].isRightExtreme())
		{
			overlapped = reverse(overlapped);
			pRightTrie->insert(overlapped, lp);
		}
		else
		{
			assert(iter->ovr.match.coord[0].isLeftExtreme());
			pLeftTrie->insert(overlapped, lp);
		}
	}		
}


// Get the complete set of overlaps for the given vertex
void SGAlgorithms::findOverlapSet(const Vertex* pVertex, VertexOverlapSet& VOSet)
{
	EdgePtrVec edges = pVertex->getEdges();

	// Add the primary overlaps to the map, and all the nodes reachable from the primaries
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pEdge = edges[i];
		Vertex* pEnd = pEdge->getEnd();
		EdgeDir transDir = pEdge->getTransitiveDir();
		Overlap ovr = pEdge->getOverlap();
		VertexOverlapPair pair = {pEnd, ovr};
		VOSet.insert(pair);

		// Recursively add nodes attached to pEnd
		SGAlgorithms::_discoverOverlaps(pVertex, pEnd, transDir, ovr, VOSet);
	}
}

// Find overlaps to vertex X via the edges of vertex Y
void SGAlgorithms::_discoverOverlaps(const Vertex* pX, const Vertex* pY, EdgeDir dir, 
                                    const Overlap& ovrXY, VertexOverlapSet& outSet)
{
	EdgePtrVec edges = pY->getEdges(dir);

	// 
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pYZ = edges[i];
		Overlap ovrYZ = pYZ->getOverlap();
		Vertex* pZ = pYZ->getEnd();
		EdgeDir dirZ = pYZ->getTransitiveDir();
	
		// Check that this vertex actually overlaps pX
		if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
		{
			Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
			VertexOverlapPair pair = { pZ, ovrXZ };
			std::pair<VertexOverlapSet::iterator, bool> ret = outSet.insert(pair);

			if(ret.second)
			{
				// The pair was actually inserted, recursively add neighbors
				_discoverOverlaps(pX, pZ, dirZ, ovrXZ, outSet);
			}
		}
	}
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
// SGRemodelVisitor - Remodel the graph to infer missing edges or remove erroneous edges
//
void SGRemodelVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

bool SGRemodelVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	bool graph_changed = false;
	(void)pGraph;
	for(size_t idx = 0; idx < ED_COUNT; idx++)
	{
		EdgeDir dir = EDGE_DIRECTIONS[idx];
		if(pVertex->countEdges(dir) > 1)
		{
			MultiOverlap mo = pVertex->getMultiOverlap();
			std::cout << "Primary MO: \n";
			mo.print();

			std::cout << "\nPrimary masked\n";
			mo.printMasked();
			
			EdgePtrVec edges = pVertex->getEdges(dir);
			for(size_t i = 0; i < edges.size(); ++i)
			{
				Edge* pXY = edges[i];
				Vertex* pY = pXY->getEnd();
				EdgeDir forwardDir = pXY->getTransitiveDir();
				EdgeDir backDir = !forwardDir;

				EdgePtrVec y_fwd_edges = pY->getEdges(forwardDir);
				EdgePtrVec y_back_edges = pY->getEdges(backDir);
				std::cout << pY->getID() << " forward edges: ";
				for(size_t j = 0; j < y_fwd_edges.size(); ++j)
					std::cout << y_fwd_edges[j]->getEndID() << ",";
				std::cout << "\n";
				
				std::cout << pY->getID() << " back edges: ";
				for(size_t j = 0; j < y_back_edges.size(); ++j)
					std::cout << y_back_edges[j]->getEndID() << ",";
				std::cout << "\n";
				std::cout << pY->getID() << " label " << pXY->getLabel() << "\n";
			}

			MultiOverlap extendedMO = SGAlgorithms::makeExtendedMultiOverlap(pVertex);
			std::cout << "\nExtended MO: \n";
			extendedMO.printMasked();

			ErrorCorrect::correctVertex(pVertex, 3, 0.01);
		}
	}
	return graph_changed;
}

//
void SGRemodelVisitor::postvisit(StringGraph*)
{
}

//
// SGErrorCorrectVisitor - Run error correction on the reads
//
bool SGErrorCorrectVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	static size_t numCorrected = 0;

	if(numCorrected > 0 && numCorrected % 50000 == 0)
		std::cerr << "Corrected " << numCorrected << " reads\n";

	std::string corrected = ErrorCorrect::correctVertex(pVertex, 5, 0.01);
	pVertex->setSeq(corrected);
	++numCorrected;
	return false;
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
	const int MIN_OVERLAP = 39;
	const double MAX_ERROR = 0.1;
	static int visited = 0;
	++visited;
	if(visited % 50000 == 0)
		std::cout << "visited: " << visited << "\n";

	// Explore the neighborhood around this graph for potentially missing overlaps
	CandidateVector candidates = getMissingCandidates(pGraph, pVertex, MIN_OVERLAP);

	MultiOverlap addedMO(pVertex->getID(), pVertex->getSeq());

	for(size_t i = 0; i < candidates.size(); ++i)
	{
		Candidate& c = candidates[i];
		int numDiff = c.ovr.match.countDifferences(pVertex->getSeq(), c.pEndpoint->getSeq());
		double error_rate = double(numDiff) / double(c.ovr.match.getMinOverlapLength());

		/*
		std::cout << "Actual:\n";
		c.ovr.match.printMatch(pVertex->getSeq(), c.pEndpoint->getSeq());
		std::cout << "ER: " << error_rate << "\n";
		*/

		if(error_rate < MAX_ERROR)
		{
			//std::cout << "OVR:\t" << c.ovr << "\n";
			Edge* p_edgeXZ = createEdges(pGraph, c.ovr, true);
			Edge* p_edgeZX = p_edgeXZ->getTwin();
			p_edgeXZ->setColor(GC_WHITE);
			p_edgeZX->setColor(GC_WHITE);
			graph_changed = true;

			addedMO.add(p_edgeXZ->getEnd()->getSeq(), p_edgeXZ->getOverlap());
		}
		/*else
		{
			std::cout << "Rejected actual:\n";
			c.ovr.match.printMatch(pVertex->getSeq(), c.pEndpoint->getSeq());
			std::cout << "ER: " << error_rate << "\n";
			std::cout << "OVR:\t" << c.ovr << "\n";
		}*/
	}

	/*
	if(graph_changed)
	{
		std::cout << "NEW MO:\n";
		addedMO.print();
	}
	*/

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
				Overlap ovrXY = pXY->getOverlap();
				Overlap ovrYZ = pYZ->getOverlap();

				if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
				{
					Overlap ovr_xz = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
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
void SGDuplicateVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
}

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

