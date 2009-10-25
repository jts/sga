//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps. See Myers (2005).
//
#include "StringGraph.h"

#define SE_CAST(x) static_cast<StringEdge*>((x))
#define CSE_CAST(x)  static_cast<const StringEdge*>((x))
#define SV_CAST(x) static_cast<StringVertex*>((x))
#define CSV_CAST(x)  static_cast<const StringVertex*>((x))

// Flip the string of an edge
void StringEdge::flip()
{
	Edge::flip();
}

// Get the edge's label (for writeDot)
std::string StringEdge::getLabel() const
{
	StringEdge* pTwin = CSE_CAST(getTwin());
	StringVertex* pEndpoint = CSV_CAST(m_pEnd);
	
	// get the unmatched coordinates in V2
	SeqCoord unmatched = pTwin->getMatchCoord().complement();
	std::string seq = unmatched.getSubstring(pEndpoint->getSeq());

	if(getComp() == EC_REVERSE)
		seq = reverseComplement(seq);
	return seq;
}

// Get the matching portion of V1 described by this edge
std::string StringEdge::getMatchStr() const
{
	StringVertex* pVert = SV_CAST(getStart());
	return m_matchCoord.getSubstring(pVert->getSeq());
}

// Return the length of the sequence
size_t StringEdge::getSeqLen() const
{
	StringEdge* pTwin = CSE_CAST(getTwin());
	SeqCoord unmatched = pTwin->getMatchCoord().complement();
	return unmatched.length();
}

// return the mapping from V1 to V2 via this edge
// If necessary, flip the SC of V2 so that it is
// in the same coordinate system as V1
Matching StringEdge::getMatch() const
{
	const StringEdge* pTwin = CSE_CAST(getTwin());

	const SeqCoord& sc = getMatchCoord();
	const SeqCoord& tsc = pTwin->getMatchCoord();

	return Matching(sc, tsc);
}

void StringEdge::validate() const
{
	/*
	const StringEdge* pTwin = CSE_CAST(getTwin());
	std::string m_v1 = getMatchStr();
	std::string m_v2 = pTwin->getMatchStr();

	if(getComp() == EC_REVERSE)
		m_v2 = reverseComplement(m_v2);

	if(m_v1 != m_v2)
	{
		std::cout << "V1M: " << m_v1 << "\n";
		std::cout << "V2M: " << m_v2 << "\n";
		std::cout << "V1MC: " << getMatchCoord() << "\n";
		std::cout << "V2MC: " << pTwin->getMatchCoord() << "\n";
		std::cout << "V1: " << SV_CAST(getStart())->getSeq() << "\n";
		std::cout << "Validation failed for edge " << *this << "\n";
		assert(false);
	}
	*/
}

// Join pEdge into the start of this edge
// This merges the intervals
void StringEdge::join(const Edge* pEdge)
{
	const StringEdge* pSE = CSE_CAST(pEdge);
	//std::cout << "J Edge : " << *pEdge << "\n";
	//std::cout << "J Edge2: " << *this << "\n";
	
	Matching m12 = pSE->getMatch();
	Matching m23 = getMatch();

	// If the V1->V2 edge is reversed, flip m12.2 to be in the same coord system as
	// V1
	if(pEdge->getComp() == EC_REVERSE)
	{
		m12.coord[1].reflect();
		m23.coord[0].reflect();
	}

	m_matchCoord = m12.inverseTranslate(m23.coord[0]);

	/*
	Matching joined;
	joined.coord[0] = m12.inverseTranslate(m23.coord[0]);
	joined.coord[1] = m23.coord[1];

	std::cout << "M12: " << m12 << "\n";
	std::cout << "M23: " << m23 << "\n";
	std::cout << "JND: " << joined << "\n";
	
	m_matchCoord = joined.coord[0];
	*/
	// Update the base class members
	Edge::join(pEdge);
}

// Join pEdge into the start of this edge
// This involves merging the intervals
void StringEdge::extend(const Edge* pEdge)
{
	// Update the base class members
	Edge::extend(pEdge);
}

void StringEdge::offsetMatch(int offset)
{
	m_matchCoord.interval.start += offset;
	m_matchCoord.interval.end += offset;
}

void StringEdge::extendMatch(int ext_len)
{
	m_matchCoord.interval.end += ext_len;
}

void StringEdge::completeMatch()
{
	if(m_matchCoord.isLeftExtreme())
		m_matchCoord.interval.end = m_matchCoord.seqlen - 1;
	else
		m_matchCoord.interval.start = 0;
}

void StringEdge::updateSeqLen(int newLen)
{
	m_matchCoord.seqlen = newLen;
}

// Merging two string vertices has two parts
// First, the sequence of the vertex is extended
// by the the content of the edge label
// Then, all the edges that are pointing to this node
// must be updated to contain the extension of the vertex
void StringVertex::merge(Edge* pEdge)
{
	// Call baseclass merge
	Vertex::merge(pEdge);

	StringEdge* pSE = SE_CAST(pEdge);
	StringEdge* pTwinSE = SE_CAST(pSE->getTwin());

	//std::cout << "Adding label to " << getID() << " str: " << pSE->getLabel() << "\n";

	// Merge the sequence
	std::string label = pSE->getLabel();
	pSE->updateSeqLen(m_seq.length() + label.length());

	bool prepend = false;

	if(pSE->getDir() == ED_SENSE)
	{
		m_seq.append(label);
	}
	else
	{
		m_seq.insert(0, label);
		prepend = true;
	}

	pSE->extendMatch(label.length());
	pTwinSE->completeMatch();

	// All the SeqCoords for the edges must have their seqlen field updated
	// Also, if we prepended sequence to this edge, all the matches in the 
	// SENSE direction must have their coordinates offset
	size_t newLen = m_seq.length();
	for(EdgePtrListIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		StringEdge* pSE2 = SE_CAST(*iter);
		pSE2->updateSeqLen(newLen);
		if(prepend && pSE2->getDir() == ED_SENSE && pSE != pSE2)
			pSE2->offsetMatch(label.length());
	}

	// Update the read count
	StringVertex* pV2 = CSV_CAST(pSE->getEnd());
	m_readCount += pV2->getReadCount();
	//validate();
}

// Ensure that the edges of the graph are correct
void StringVertex::validate() const
{
	Vertex::validate();

	for(EdgePtrListConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		StringEdge* pSE = SE_CAST(*iter);
		pSE->validate();
		/*
		std::string label = pSE->getLabel();
		StringVertex* pEnd = SV_CAST(pSE->getEnd());

		std::string vertSeq = pEnd->getSeq();
		EdgeDir suffixDir = (pSE->getComp() == EC_SAME) ? pSE->getDir() : !pSE->getDir();
		std::string vertSuffix = (suffixDir == ED_SENSE) ? vertSeq.substr(vertSeq.length() - label.length()) :
																vertSeq.substr(0, label.length());

		if(pSE->getComp() == EC_REVERSE)
			vertSuffix = reverseComplement(vertSuffix);
		if(vertSuffix != label)
		{
			std::cerr << "Warning edge label " << label << " does not match vertex suffix " << vertSuffix << "\n";
		}
		*/
	}
}

void StringVertex::sortAdjList()
{
	StringEdgeLenComp comp;
	m_edges.sort(comp);
}

// Compare string edge points by length
bool StringEdgeLenComp::operator()(const Edge* pA, const Edge* pB)
{
	const StringEdge* pSA = static_cast<const StringEdge*>(pA);
	const StringEdge* pSB = static_cast<const StringEdge*>(pB);
	return pSA->getSeqLen() < pSB->getSeqLen();
}

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
		EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted

		if(edges.size() == 2)
		{
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
						++num_bubbles;
						bubble_found = true;
					}
					else
					{
						// Endpoint has not been hit, set it to visited
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



