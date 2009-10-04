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

// Flip the string of an edge
void StringEdge::flip()
{
	Edge::flip();
	m_seq = reverseComplement(m_seq);
}

// Get the edge's label (for writeDot)
std::string StringEdge::getLabel() const
{
	return m_seq;
}

// StringEdge::merge is inherited from Edge::merge
// The baseclass merge updates the endpoint of the
// Edge, this function updates the data
// which is the concatenation of the edge sequences
void StringEdge::merge(const Edge* pEdge)
{
	// Call baseclass
	Edge::merge(pEdge);

	// Update sequence
	const StringEdge* pSE = static_cast<const StringEdge*>(pEdge);
	updateLabel(pSE);
}

// Update the label for this edge to include the label of the provided edge
void StringEdge::updateLabel(const StringEdge* pSE)
{
	std::string edgeLabel = pSE->getSeq();
	// The edge label must be flipped if the edge from pSE is
	// not in the same direction as this edge. This implies that 
	// the orientation of V2 and V3 (the endpoints of pSE) are different orientation
	if(getDir() != pSE->getDir())
		edgeLabel = reverseComplement(edgeLabel);
	if(getDir() == ED_SENSE)
		m_seq.append(edgeLabel);
	else
		m_seq.insert(0, edgeLabel);
	//std::cout << "Updated edge " << *this << " to be: " << m_seq << "\n";
}

// Merging two string vertices has two parts
// First, the sequence of the vertex is extended
// by the the content of the edge label
// Then, all the edges that are pointing to this node
// must be updated to contain the extension of the vertex
void StringVertex::merge(const Edge* pEdge)
{
	// Call baseclass merge
	Vertex::merge(pEdge);

	const StringEdge* pSE = static_cast<const StringEdge*>(pEdge);
	StringVertex* pV2 = static_cast<StringVertex*>(pSE->getEnd());

	// Merge the sequence
	if(pSE->getDir() == ED_SENSE)
		m_seq.append(pSE->getSeq());
	else
		m_seq.insert(0, pSE->getSeq()); //prepend

	// Now update the edges that point TO this vertex to contain the merged-in label
	for(EdgePtrListIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		// Only update edges in the opposite direction of the merged-in edge
		// This is the dimension that grew
		if((*iter)->getDir() != pEdge->getDir())
		{
			StringEdge* pTwinSE = static_cast<StringEdge*>((*iter)->getTwin());
			pTwinSE->updateLabel(pSE);
		}
	}

	// Update the read count
	m_readCount += pV2->getReadCount();
	//validate();
}

// Ensure that the edges of the graph are correct
void StringVertex::validate() const
{
	Vertex::validate();

	for(EdgePtrListConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		StringEdge* pSE = static_cast<StringEdge*>(*iter);
		std::string label = pSE->getSeq();
		StringVertex* pEnd = static_cast<StringVertex*>(pSE->getEnd());

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
}

// Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
bool SGTransRedVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
	(void)pGraph;
	(void)pVertex;
	
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
		size_t longestLen = pLongestEdge->getSeqLen();
		
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
							//std::cout << "Marking " << pWXEdge->getEndID() << " as transitive\n";
						}
					}
					else
						break;
				}
			}
		}

		for(size_t i = 0; i < edges.size(); ++i)
		{
			if(edges[i]->getEnd()->getColor() == GC_BLACK)
			{
				// Mark the edge and its twin for removal
				edges[i]->setColor(GC_BLACK);
				edges[i]->getTwin()->setColor(GC_BLACK);
			}
			edges[i]->getEnd()->setColor(GC_WHITE);
		}
	}

	return false;
}

// Remove all the marked edges
void SGTransRedVisitor::postvisit(StringGraph* pGraph)
{
	pGraph->sweepEdges(GC_BLACK);
}

