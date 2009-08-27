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
	for(EdgePtrMapIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		// Only update edges in the opposite direction of the merged-in edge
		// This is the dimension that grew
		if(iter->second->getDir() != pEdge->getDir())
		{
			StringEdge* pTwinSE = static_cast<StringEdge*>(iter->second->getTwin());
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

	for(EdgePtrMapConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		StringEdge* pSE = static_cast<StringEdge*>(iter->second);
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

// Visitor which outputs the graph in fasta format
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	StringVertex* pSV = static_cast<StringVertex*>(pVertex);
	m_fileHandle << ">" << pSV->getID() << " " <<  pSV->getReadCount() << "\n";
	m_fileHandle << pSV->getSeq() << "\n";
	return false;
}

