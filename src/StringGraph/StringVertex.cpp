//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// StringVertex - Derived from Bigraph/Vertex
//
#include "StringVertex.h"
#include "StringEdge.h"

// 
StringVertex::~StringVertex()
{
	// Inform the pair of this vertex that it is being deleted
	if(m_pPairVertex != NULL)
		m_pPairVertex->clearPairVertex();
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
	size_t label_len = label.size();
	pSE->updateSeqLen(m_seq.length() + label.length());
	bool prepend = false;

	if(pSE->getDir() == ED_SENSE)
	{
		m_seq.append(label);
	}
	else
	{
		label.insert(label.size(), m_seq);
		std::swap(m_seq, label);
		prepend = true;
	}

	pSE->extendMatch(label_len);
	pTwinSE->extendMatchFullLength();

	// All the SeqCoords for the edges must have their seqlen field updated
	// Also, if we prepended sequence to this edge, all the matches in the 
	// SENSE direction must have their coordinates offset
	size_t newLen = m_seq.length();
	for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		StringEdge* pSE2 = SE_CAST(*iter);
		pSE2->updateSeqLen(newLen);
		if(prepend && pSE2->getDir() == ED_SENSE && pSE != pSE2)
			pSE2->offsetMatch(label_len);
	}

	// Update the read count
	const StringVertex* pV2 = CSV_CAST(pSE->getEnd());
	m_readCount += pV2->getReadCount();

#ifdef VALIDATE
	VALIDATION_WARNING("StringVertex::merge")
	validate();
#endif

}

// Ensure that the edges of the graph are correct
void StringVertex::validate() const
{
	Vertex::validate();

	for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
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
	std::sort(m_edges.begin(), m_edges.end(), comp);
}

// Ensure that this vertex has at most one edge per direction to any other vertex
void StringVertex::makeUnique()
{
	// Sort the edge lists by length
	sortAdjList();
	EdgePtrVec uniqueVec;
	makeUnique(ED_SENSE, uniqueVec);
	makeUnique(ED_ANTISENSE, uniqueVec);
	m_edges.swap(uniqueVec);
}

// 
void StringVertex::makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec)
{
	std::set<VertexID> idSet;
	for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		Edge* pEdge = *iter;
		if(pEdge->getDir() == dir)
		{
			std::pair<std::set<VertexID>::iterator, bool> result = idSet.insert(pEdge->getEndID());
			if(result.second == true)
			{
				uniqueVec.push_back(*iter);
			}
			else
			{
				std::cerr << getID() << " has a duplicate edge to " << pEdge->getEndID() << " in direction " << dir << "\n";
				
				// Delete the edge and remove it from the twin
				Edge* pTwin = pEdge->getTwin();
				Vertex* pPartner = pEdge->getEnd();
				pPartner->removeEdge(pTwin);
				delete pEdge;
				pEdge = NULL;
				delete pTwin;
				pTwin = NULL;
			}
		}
	}
}

