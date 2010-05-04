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

// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
// The algorithm is as follows. All the vertices that can be reached by vertices other
// than the deletion vertex are marked as reachable and not considered. The vertices
// uniquely reachable through the deletion vertex are considered in order of putative
// overlap with pVertex. If a new edge is created, all vertices reachable from that
// edge are subsequently marked as reachable so that no transitive edges are created.
void SGAlgorithms::remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge)
{
	assert(pVertex == pDeleteEdge->getStart());
	EdgePtrVec edges = pVertex->getEdges();

	// this is the initial overlap between pVertex and pDeletionVertex
	Overlap ovrXY = pDeleteEdge->getOverlap();
	EdgeDesc edXY = pDeleteEdge->getDesc();

	// Construct the set of vertices that are reachable by valid edges
	EdgeDescOverlapMap exclusionSet;
	exclusionSet.insert(std::make_pair(edXY, ovrXY));

	// Recursively add the vertices connected to pX to the exclusionSet
	// except for the neighbors that are exclusively reachable from pDeleteVertex
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pEdge = edges[i];
		
		// Skip adding vertices for the deletion edge
		if(pEdge != pDeleteEdge)
		{
			Overlap ovr = pEdge->getOverlap();
			EdgeDesc ed = pEdge->getDesc();
			exclusionSet.insert(std::make_pair(ed, ovr));

			// Recursively add the neighbors of pEnd to the set
			addOverlapsToSet(pVertex, ed, ovr, exclusionSet);
		}
	}
	
	// Build the initial set of potential new overlaps from the 
	// neighbors of pDeleteVertex.
	EdgeDescSet seenEdges;
	ExplorePriorityQueue exploreQueue;
	enqueueEdges(pVertex, edXY, ovrXY, exploreQueue, seenEdges, &exclusionSet);

	// Iterate through the queue in order of overlap length
	while(!exploreQueue.empty())
	{
		ExploreElement currElem = exploreQueue.top();
		exploreQueue.pop();
		std::cout << "Current element: " << currElem.ovr << " length: " << currElem.ovr.match.coord[0].length() << "\n";
		// Case 1, endpoint is reachable from some other edge of pVertex
		// and is therefore transitive
		if(exclusionSet.count(currElem.ed) > 0)
		{
			std::cout << "	Elem is in reachable, skipping\n";
			continue;
		}

		// Case 2, this may form a valid edge
		double error_rate = calcErrorRate(pVertex, currElem.ed.pVertex, currElem.ovr);
		int overlap_len = currElem.ovr.match.getMinOverlapLength();
		if(overlap_len >= pGraph->getMinOverlap())
		{
			WARN_ONCE("CHECK ERROR RATE WITHIN EPSILON")
			if(error_rate - 0.0001 <= pGraph->getErrorRate())
			{
				std::cout << "Adding edge with overlap " << currElem.ovr << "\n";
				Edge* pCreatedEdge = SGUtil::createEdges(pGraph, currElem.ovr, false);
				std::cout << "CreatedED: " << pCreatedEdge->getDesc() << " expect: " << currElem.ed << "\n";
				assert(pCreatedEdge != NULL);
				assert(pCreatedEdge->getDesc() == currElem.ed);
				
				// This vertex is now connected to pVertex, add its neighbors to the exclusion set
				addOverlapsToSet(pVertex, currElem.ed, currElem.ovr, exclusionSet);
			}
			else
			{
				std::cout << "	Error rate too high: " << error_rate << "\n";
			}
		}
	}
}

// Add the neighbors of pY to the explore queue if they overlap pX. If pExclusionSet
// is not NULL, do not add neighbors of pY to the queue
void SGAlgorithms::enqueueEdges(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
                                ExplorePriorityQueue& outQueue, EdgeDescSet& seenEdges, 
								EdgeDescOverlapMap* pExclusionSet)
{
	Vertex* pY = edXY.pVertex;
	EdgeDir dirY = correctDir(edXY.dir, edXY.comp);
	EdgePtrVec neighborEdges = pY->getEdges(dirY);

	for(size_t i = 0; i < neighborEdges.size(); ++i)
	{
		Edge* pEdgeYZ = neighborEdges[i];
		if(pEdgeYZ->getEnd() != pX)
		{
			EdgeDesc edYZ = pEdgeYZ->getDesc();
			EdgeDesc edXZ = SGAlgorithms::inferTransitiveEdgeDesc(edXY, edYZ);
			bool isExcluded = pExclusionSet != NULL && pExclusionSet->count(edXZ) > 0;
			if(!isExcluded && seenEdges.count(edXZ) == 0)
			{
				Overlap ovrYZ = pEdgeYZ->getOverlap();
				
				// Check that this vertex actually overlaps pX
				if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
				{
					Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
					std::cout << "Inferred overlap: " << ovrXZ << " ed: " << edXZ << " from: " << pY->getID() << "\n";
					ExploreElement elem(edXZ, ovrXZ);
					outQueue.push(elem);
					seenEdges.insert(edXZ);
					enqueueEdges(pX, edXZ, ovrXZ, outQueue, seenEdges, pExclusionSet);
				}
			}
			else
			{
				/*
				Overlap ovrYZ = pEdgeYZ->getOverlap();
				Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
				std::cout << "Rejected enqueue: " << ovrXZ << "\n";
				*/
			}
		}
	}
}

// Recursively add overlaps to pX inferred from the edges of pY to outMap
void SGAlgorithms::addOverlapsToSet(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, EdgeDescOverlapMap& outMap)
{
	WARN_ONCE("HANDLE CONTAINMENT EDGES?");

	Vertex* pY = edXY.pVertex;

	// Calculate the direction of the edge coming out of Y
	EdgeDir dirY = correctDir(edXY.dir, edXY.comp);
	EdgePtrVec neighborEdges = pY->getEdges(dirY);
	for(size_t i = 0; i < neighborEdges.size(); ++i)
	{
		Edge* pEdgeYZ = neighborEdges[i];
		EdgeDesc edYZ = pEdgeYZ->getDesc();
		EdgeDesc edXZ = SGAlgorithms::inferTransitiveEdgeDesc(edXY, edYZ);
		Vertex* pZ = pEdgeYZ->getEnd();

		if(pZ != pX && outMap.count(edXZ) == 0)
		{
			Overlap ovrYZ = pEdgeYZ->getOverlap();
			EdgeDesc edYZ = pEdgeYZ->getDesc();

			// Check that this vertex actually overlaps pX
			if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
			{
				Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
				EdgeDesc edXZ = SGAlgorithms::inferTransitiveEdgeDesc(edXY, edYZ);

				outMap.insert(std::make_pair(edXZ, ovrXZ));
				std::cout << "Vertex " << pZ->getID() << " is reachable from " << pY->getID() << " " << ovrXZ << "\n";
				std::cout << "EdgeDesc XY: " << edXY << " edYZ: " << edYZ << " edXZ: " << edXZ << "\n";
				addOverlapsToSet(pX, edXZ, ovrXZ, outMap);
			}
		}
	}
}

// Calculate the error rate between the two vertices
double SGAlgorithms::calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY)
{
	int num_diffs = ovrXY.match.countDifferences(pX->getSeq(), pY->getSeq());
	return static_cast<double>(num_diffs) / static_cast<double>(ovrXY.match.getMinOverlapLength());
}

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

// Infer an EdgeDesc between X -> Z
// given EdgeDescs X -> Y, Y -> Z
// The input edges are between X->Y Y->Z
// and the returned overlap is X->Z
EdgeDesc SGAlgorithms::inferTransitiveEdgeDesc(const EdgeDesc& edXY, const EdgeDesc& edYZ)
{
	EdgeDesc out;
	out.pVertex = edYZ.pVertex; // the endpoint is Z
	out.dir = edXY.dir; // it must be in the same direction as X->Y
	out.comp = (edYZ.comp == EC_REVERSE) ? !edXY.comp : edXY.comp;
	return out;
}

// Return true if XZ has an overlap
bool SGAlgorithms::hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
	Match match_yx = ovrXY.match;
	match_yx.swap(); 
	Match match_yz = ovrYZ.match;
	return Match::doMatchesIntersect(match_yx, match_yz);
}

// Construct an extended multioverlap for a vertex
MultiOverlap SGAlgorithms::makeExtendedMultiOverlap(const Vertex* pVertex)
{
	EdgeDescOverlapMap overlapMap;
	findOverlapMap(pVertex, overlapMap);

	MultiOverlap mo(pVertex->getID(), pVertex->getSeq());
	for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin();
	    iter != overlapMap.end(); ++iter)
	{
		mo.add(iter->first.pVertex->getSeq(), iter->second);
	}
	return mo;
}

//
void SGAlgorithms::makeExtendedSeqTries(const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie)
{
	double lp = log(p_error);
	EdgeDescOverlapMap overlapMap;
	findOverlapMap(pVertex, overlapMap);

	for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin();
	    iter != overlapMap.end(); ++iter)
	{
		// Coord[0] of the match is wrt pVertex, coord[1] is the other read
		std::string overlapped = iter->second.match.coord[1].getSubstring(iter->first.pVertex->getSeq());
		if(iter->second.match.isRC())
			overlapped = reverseComplement(overlapped);

		if(iter->second.match.coord[0].isRightExtreme())
		{
			overlapped = reverse(overlapped);
			pRightTrie->insert(overlapped, lp);
		}
		else
		{
			assert(iter->second.match.coord[0].isLeftExtreme());
			pLeftTrie->insert(overlapped, lp);
		}
	}		
}


// Get the complete set of overlaps for the given vertex
void SGAlgorithms::findOverlapMap(const Vertex* pVertex, EdgeDescOverlapMap& outMap)
{
	EdgePtrVec edges = pVertex->getEdges();

	// Add the primary overlaps to the map, and all the nodes reachable from the primaries
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pEdge = edges[i];
		EdgeDesc ed = pEdge->getDesc();
		Overlap ovr = pEdge->getOverlap();
		outMap.insert(std::make_pair(ed, ovr));
		std::cout << "Vertex " << ed.pVertex->getID() << " is reachable from " << pVertex->getID() << "\n";

		// Recursively add nodes attached to pEnd
		SGAlgorithms::_discoverOverlaps(pVertex, ed, ovr, outMap);
	}
}

// Find overlaps to vertex X via the edges of vertex Y
void SGAlgorithms::_discoverOverlaps(const Vertex* pX, const EdgeDesc& edXY, 
                                    const Overlap& ovrXY, EdgeDescOverlapMap& outMap)
{
	Vertex* pY = edXY.pVertex;
	EdgeDir dirY = correctDir(edXY.dir, edXY.comp);
	EdgePtrVec edges = pY->getEdges(dirY);

	// 
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Edge* pYZ = edges[i];
		Overlap ovrYZ = pYZ->getOverlap();
		EdgeDesc edYZ = pYZ->getDesc();
		Vertex* pZ = pYZ->getEnd();

		// Ignore self overlaps
		if(pZ == pX)
			continue;

		// Check that this vertex actually overlaps pX
		if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
		{
			Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
			EdgeDesc edXZ = SGAlgorithms::inferTransitiveEdgeDesc(edXY, edYZ);

			std::pair<EdgeDescOverlapMap::iterator, bool> ret = outMap.insert(std::make_pair(edXZ, ovrXZ));
			if(ret.second)
			{
				std::cout << "Vertex " << edXZ.pVertex->getID() << " is reachable from " << pY->getID() << "\n";
				// The pair was inserted, recursively add neighbors
				_discoverOverlaps(pX, edXZ, ovrXZ, outMap);
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
	size_t trans_count = 0;
	static const size_t FUZZ = 200; // see myers...

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
						//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << " in stage 2\n";
					}
				}
				else
					break;
			}
		}

		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_BLACK)
			{
				// Mark the edge and its twin for removal
				if(edges[i]->getColor() != GC_BLACK || edges[i]->getTwin()->getColor() != GC_BLACK)
				{
					edges[i]->setColor(GC_BLACK);
					edges[i]->getTwin()->setColor(GC_BLACK);
					marked_edges += 2;
					trans_count++;
				}
			}
			edges[i]->getEnd()->setColor(GC_WHITE);
		}
	}

	if(trans_count > 0)
		++marked_verts;

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
bool SGContainRemoveVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	// Skip the computation if this vertex has already been marked
	if(pVertex->getColor() == GC_BLACK)
		return false;

	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Overlap ovr = edges[i]->getOverlap();
		Match m = edges[i]->getMatch();
		if(ovr.match.isContainment())
		{
			Vertex* pVertexY = edges[i]->getEnd();
			// Skip the resolution step if the vertex has already been marked
			if(pVertexY->getColor() == GC_BLACK)
				continue;

			Vertex* pToRemove = NULL;
			
			// If containedIdx is 0, then this vertex is the one to remove
			if(ovr.getContainedIdx() == 0)
			{
				pToRemove = pVertex;
			}
			else
			{
				pToRemove = pVertexY;
			}

			assert(pToRemove != NULL);
			pToRemove->setColor(GC_BLACK);

			std::cout << "REMOVING " << pToRemove->getID() << "\n";
			// Add any new irreducible edges that exist when pToRemove is deleted
			// from the graph
			EdgePtrVec neighborEdges = pToRemove->getEdges();
			for(size_t j = 0; j < neighborEdges.size(); ++j)
			{
				Vertex* pRemodelVert = neighborEdges[j]->getEnd();
				Edge* pRemodelEdge = neighborEdges[j]->getTwin();
				std::cout << "Patching " << pRemodelVert->getID() << "\n";
				SGAlgorithms::remodelVertexForExcision(pGraph, 
				                                       pRemodelVert, 
													   pRemodelEdge);
			}
			
			// Delete the edges from the graph
			for(size_t j = 0; j < neighborEdges.size(); ++j)
			{
				Vertex* pRemodelVert = neighborEdges[j]->getEnd();
				Edge* pRemodelEdge = neighborEdges[j]->getTwin();
				pRemodelVert->deleteEdge(pRemodelEdge);
				pToRemove->deleteEdge(neighborEdges[j]);
			}
			
			(void)pGraph;
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
// SGEdgeStatsVisitor - Compute and display summary statistics of 
// the overlaps in the graph, including edges that were potentially missed
//
void SGEdgeStatsVisitor::previsit(StringGraph* pGraph)
{
	pGraph->setColors(GC_WHITE);
	maxDiff = 0;
	minOverlap = pGraph->getMinOverlap();
	maxOverlap = 0;

}

bool SGEdgeStatsVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	const int MIN_OVERLAP = pGraph->getMinOverlap();
	const double MAX_ERROR = pGraph->getErrorRate();

	static int visited = 0;
	++visited;
	if(visited % 50000 == 0)
		std::cout << "visited: " << visited << "\n";

	// Add stats for the found overlaps
	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Overlap ovr = edges[i]->getOverlap();
		int numDiff = ovr.match.countDifferences(pVertex->getSeq(), edges[i]->getEnd()->getSeq());
		int overlapLen = ovr.match.getMinOverlapLength();
		addOverlapToCount(overlapLen, numDiff, foundCounts);
	}
		
	// Explore the neighborhood around this graph for potentially missing overlaps
	CandidateVector candidates = getMissingCandidates(pGraph, pVertex, MIN_OVERLAP);
	MultiOverlap addedMO(pVertex->getID(), pVertex->getSeq());
	for(size_t i = 0; i < candidates.size(); ++i)
	{
		Candidate& c = candidates[i];
		int numDiff = c.ovr.match.countDifferences(pVertex->getSeq(), c.pEndpoint->getSeq());
		double error_rate = double(numDiff) / double(c.ovr.match.getMinOverlapLength());

		if(error_rate < MAX_ERROR)
		{
			int overlapLen = c.ovr.match.getMinOverlapLength();
			addOverlapToCount(overlapLen, numDiff, missingCounts);
		}
	}
	
	return false;
}

//
void SGEdgeStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{	
	printf("FoundOverlaps\n");
	printCounts(foundCounts);

	printf("\nPotentially Missing Overlaps\n\n");
	printCounts(missingCounts);
}

//
void SGEdgeStatsVisitor::printCounts(CountMatrix& matrix)
{
	// Header row
	printf("OL\t");
	for(int j = 0; j <= maxDiff; ++j)
	{
		printf("%d\t", j);
	}

	printf("sum\n");
	IntIntMap columnTotal;
	for(int i = minOverlap; i <= maxOverlap; ++i)
	{
		printf("%d\t", i);
		int sum = 0;
		for(int j = 0; j <= maxDiff; ++j)
		{
			int v = matrix[i][j];
			printf("%d\t", v);
			sum += v;
			columnTotal[j] += v;
		}
		printf("%d\n", sum);
	}

	printf("total\t");
	int total = 0;
	for(int j = 0; j <= maxDiff; ++j)
	{
		int v = columnTotal[j];
		printf("%d\t", v);
		total += v;
	}
	printf("%d\n", total);
}

//
void SGEdgeStatsVisitor::addOverlapToCount(int ol, int nd, CountMatrix& matrix)
{
	matrix[ol][nd]++;

	if(nd > maxDiff)
		maxDiff = nd;

	if(ol > maxOverlap)
		maxOverlap = ol;
}

// Explore the neighborhood around a vertex looking for missing overlaps
SGEdgeStatsVisitor::CandidateVector SGEdgeStatsVisitor::getMissingCandidates(StringGraph* /*pGraph*/, 
                                                                             Vertex* pVertex, 
																			 int minOverlap) const
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
	sum_edgeLen = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGGraphStatsVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	int s_count = pVertex->countEdges(ED_SENSE);
	int as_count = pVertex->countEdges(ED_ANTISENSE);
	if(s_count == 0 && as_count == 0)
		++num_island;
	else if(s_count == 0 || as_count == 0)
	{
		std::cout << "TERMINAL: " << pVertex->getID() << "\n";
		++num_terminal;
	}

	if(s_count > 1 && as_count > 1)
		++num_dibranch;
	else if(s_count > 1 || as_count > 1)
		++num_monobranch;

	if(s_count == 1 || as_count == 1)
		++num_transitive;

	num_edges += (s_count + as_count);
	++num_vertex;

	EdgePtrVec edges = pVertex->getEdges();
	for(size_t i = 0; i < edges.size(); ++i)
		sum_edgeLen += edges[i]->getSeqLen();

	return false;
}

// Remove all the marked edges
void SGGraphStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{
	printf("island: %d terminal: %d monobranch: %d dibranch: %d transitive: %d\n", num_island, num_terminal,
	                                                                               num_monobranch, num_dibranch, num_transitive);
	printf("Total Vertices: %d Total Edges: %d Sum edge length: %zu\n", num_vertex, num_edges, sum_edgeLen);
}

