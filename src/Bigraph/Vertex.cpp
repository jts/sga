//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#include "Vertex.h"
#include "Edge.h"
#include <algorithm>

Vertex::~Vertex()
{
	EdgePtrVecIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		delete *iter;
		*iter = NULL;
	}

	if(m_pPairVertex != NULL)
		m_pPairVertex->clearPairVertex();
}

// Merging two string vertices has two parts
// First, the sequence of the vertex is extended
// by the the content of the edge label
// Then, all the edges that are pointing to this node
// must be updated to contain the extension of the vertex
void Vertex::merge(Edge* pEdge)
{
	Edge* pTwin = pEdge->getTwin();
	//std::cout << "Adding label to " << getID() << " str: " << pSE->getLabel() << "\n";

	// Merge the sequence
	std::string label = pEdge->getLabel();
	size_t label_len = label.size();
	pEdge->updateSeqLen(m_seq.length() + label.length());
	bool prepend = false;

	if(pEdge->getDir() == ED_SENSE)
	{
		m_seq.append(label);
	}
	else
	{
		label.insert(label.size(), m_seq);
		std::swap(m_seq, label);
		prepend = true;
	}

	pEdge->extendMatch(label_len);
	pTwin->extendMatchFullLength();

	// All the SeqCoords for the edges must have their seqlen field updated
	// Also, if we prepended sequence to this edge, all the matches in the 
	// SENSE direction must have their coordinates offset
	size_t newLen = m_seq.length();
	for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		Edge* pUpdateEdge = *iter;
		pUpdateEdge->updateSeqLen(newLen);
		if(prepend && pUpdateEdge->getDir() == ED_SENSE && pEdge != pUpdateEdge)
			pUpdateEdge->offsetMatch(label_len);
	}

	// Update the read count
	m_readCount += pEdge->getEnd()->getReadCount();

#ifdef VALIDATE
	VALIDATION_WARNING("Vertex::merge")
	validate();
#endif

}

void Vertex::validate() const
{
	Vertex::validate();

	for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		(*iter)->validate();
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

//
void Vertex::sortAdjListByID()
{
	EdgeIDComp comp;
	std::sort(m_edges.begin(), m_edges.end(), comp);
}

void Vertex::sortAdjListByLen()
{
	EdgeLenComp comp;
	std::sort(m_edges.begin(), m_edges.end(), comp);
}

// Ensure that this vertex has at most one edge per direction to any other vertex
void Vertex::makeUnique()
{
	// Sort the edge lists by length
	sortAdjListByLen();
	EdgePtrVec uniqueVec;
	makeUnique(ED_SENSE, uniqueVec);
	makeUnique(ED_ANTISENSE, uniqueVec);
	m_edges.swap(uniqueVec);
}

// 
void Vertex::makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec)
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


// Add an edge
void Vertex::addEdge(Edge* ep)
{
	assert(ep->getStart() == this);

#ifdef VALIDATE
	for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getEndID() == ep->getEndID())
		{
			std::cout << "Attempted to add duplicate edge with ID: " << ep->getEndID() 
						<< " to vertex: " << ep->getStartID() << "\n";
			std::cout << "Added in desc: " << ep->getDesc() << " curr desc: " << (*iter)->getDesc() << "\n";
			assert(false);
		}
	}
#endif
	m_edges.push_back(ep);
}

// Remove an edge
void Vertex::removeEdge(Edge* pEdge)
{
	// Check if the edge exists
	assert(pEdge != NULL);
	removeEdge(pEdge->getDesc());
}

// Remove edge
// Note - this does not delete the edge through the pointer
void Vertex::removeEdge(const EdgeDesc& ed)
{

	EdgePtrVecIter iter = findEdge(ed);
	assert(iter != m_edges.end());
	m_edges.erase(iter);
}

// Delete all the edges, and their twins, from this vertex
void Vertex::deleteEdges()
{
	EdgePtrVecIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		Edge* pEdge = *iter;
		Edge* pTwin = pEdge->getTwin();
		Vertex* pPartner = pEdge->getEnd();
		pPartner->removeEdge(pTwin);
		delete pEdge;
		pEdge = NULL;
		delete pTwin;
		pTwin = NULL;
		*iter = NULL;
	}
	m_edges.clear();
}

// Delete edges that are marked
// This only deletes the edge and not its twin
void Vertex::sweepEdges(GraphColor c)
{
	EdgePtrVecIter iter = m_edges.begin();
	while(iter != m_edges.end())
	{
		Edge* pEdge = *iter;
		if(pEdge->getColor() == c)
		{
			delete pEdge;
			pEdge = NULL;
			iter = m_edges.erase(iter);
		}
		else
			++iter;
	}
}

// Return the iterator to the edge matching edgedesc
EdgePtrVecIter Vertex::findEdge(const EdgeDesc& ed)
{
	for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDesc() == ed)
			return iter;
	}
	return m_edges.end();
}

//
EdgePtrVecConstIter Vertex::findEdge(const EdgeDesc& ed) const
{
	for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDesc() == ed)
			return iter;
	}
	return m_edges.end();
}

// Check for the presence of an edge
bool Vertex::hasEdge(Edge* pEdge) const
{
	return hasEdge(pEdge->getDesc());
}

//
bool Vertex::hasEdge(const EdgeDesc& ed) const
{
	return findEdge(ed) != m_edges.end();
}

// Return the edge matching the descriptions
Edge* Vertex::getEdge(const EdgeDesc& ed)
{
	 EdgePtrVecIter i = findEdge(ed);
	 assert(i != m_edges.end());
	 return *i;
}

// Find edges to the specified vertex
EdgePtrVec Vertex::findEdgesTo(VertexID id)
{
	EdgePtrVecConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
	{
		if((*iter)->getEndID() == id)
			outEdges.push_back(*iter);
	}
	return outEdges;
}

//
void Vertex::setPairVertex(Vertex* pPair) 
{ 
	m_pPairVertex = pPair; 
}

//
void Vertex::clearPairVertex() 
{ 
	m_pPairVertex = NULL; 
}


//
// Get the edges in a particular direction
// This preserves the ordering of the edges
//
EdgePtrVec Vertex::getEdges(EdgeDir dir)
{
	EdgePtrVecConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDir() == dir)
			outEdges.push_back(*iter);
	}
	return outEdges;
}


// Get the edges
EdgePtrVec Vertex::getEdges()
{
	EdgePtrVecConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
		outEdges.push_back(*iter);
	return outEdges;	
}

void Vertex::setEdgeColors(GraphColor c) 
{ 
	for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
		(*iter)->setColor(c);
}

// Count the edges
// This function is not necessarily constant time
size_t Vertex::countEdges() const
{ 
	return m_edges.size(); 
}

//
size_t Vertex::countEdges(EdgeDir dir)
{
	EdgePtrVec ev = getEdges(dir);
	return ev.size();
}

// Return the amount of memory this vertex is using, in bytes
size_t Vertex::getMemSize() const
{
	return sizeof(*this) + (m_edges.size() * sizeof(Edge*));
}


// Output edges in graphviz format
void Vertex::writeEdges(std::ostream& out, int dotFlags) const
{
	EdgePtrVecConstIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		if(dotFlags & DF_UNDIRECTED)
		{
			if((*iter)->getStartID() < (*iter)->getEndID())
			{
				out << "\"" << (*iter)->getStart() << "\" -- \"" << (*iter)->getEnd() << "\"";
			}
		}
		else
		{
			out << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
			std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
			std::string label = ((*iter)->getComp() == EC_SAME) ? "S" : "F";
			out << "\" [color=\"" << color << "\" ";
			out << "label=\"" << (*iter)->getLabel() << "\"];";
		}
		out << "\n";
	}
}

bool EdgeIDComp::operator()(const Edge* pA, const Edge* pB) 
{
   	return pA->getEndID() < pB->getEndID();
}

// Compare string edge points by length
bool EdgeLenComp::operator()(const Edge* pA, const Edge* pB)
{
	return pA->getSeqLen() < pB->getSeqLen();
}

