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
#include "SeqReader.h"

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
	if(pSE->getComp() == EC_REVERSE)
		edgeLabel = reverseComplement(edgeLabel);
	if(getDir() == ED_SENSE)
		m_seq.append(edgeLabel);
	else
		m_seq.insert(0, edgeLabel);
	std::cout << "Updated edge to be: " << m_seq << "\n";
}


// StringVertex::Merge is a virtual function inheriting from Vertex::merge
// which joins two vertices together
// The sequence of pEdge is the unmatched portion of V2 that 
// needs to be added to V1
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
		// (This is the dimension of the sequence that grew, therefore the edge label
		// must grow)
		if(iter->second->getDir() != pEdge->getDir())
		{
			StringEdge* pTwinSE = static_cast<StringEdge*>(iter->second->getTwin());
			pTwinSE->updateLabel(pSE);
		}
	}

	// Update the read count
	m_readCount += pV2->getReadCount();

	// As a sanity check, ensure that v2->seq is now contained by 
	// v1->seq
	std::cout << "this:  " << m_seq << "\n"; 
	std::cout << "other: " << pV2->getSeq() << "\n";
	assert(getSeq().find(pV2->getSeq()) != std::string::npos);
}

// Visitor which outputs the graph in fasta format
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
	StringVertex* pSV = static_cast<StringVertex*>(pVertex);
	m_fileHandle << ">" << pSV->getID() << " " <<  pSV->getReadCount() << "\n";
	m_fileHandle << pSV->getSeq() << "\n";
	return false;
}

// Construct a string graph
StringGraph* createStringGraph(std::string readFile, std::string overlapFile)
{
	StringGraph* pGraph = new StringGraph;

	// Add the reads as the vertices
	SeqReader reader(readFile);
	SeqItem si;
	while(reader.get(si))
		pGraph->addVertex(new StringVertex(si.id, si.seq));

	// Add the overlaps as edges
	std::ifstream overlapReader(overlapFile.c_str());
	Overlap o;
	while(overlapReader >> o)
	{
		std::cout << "Read: " << o << "\n";

		// Initialize data and perform checks
		StringVertex* pVerts[2];
		size_t sizes[2];
		Sequence overhangs[2];

		for(size_t idx = 0; idx < 2; ++idx)
		{
			pVerts[idx] = static_cast<StringVertex*>(pGraph->getVertex(o.read[idx].id));
			assert(pVerts[idx]);
			sizes[idx] = pVerts[idx]->getSeq().length();
			
			// Ensure that the overlaps are canonical
			assert(o.read[idx].interval.start < o.read[idx].interval.end);
			
			// Ensure the reads are not identical
			assert(!o.read[idx].isContainment(sizes[idx]));

			overhangs[idx] = getOverhangString(o.read[idx], pVerts[idx]->getSeq());
		}

		// Ensure that each overlap is extreme
		if(!o.read[0].isExtreme(sizes[0]) || !o.read[1].isExtreme(sizes[1]))
			continue;

		EdgeComp comp = (o.read[0].isReverse() || o.read[1].isReverse()) ? EC_REVERSE : EC_SAME;

		// Add edges
		StringEdge* pEdges[2];
		for(size_t idx = 0; idx < 2; ++idx)
		{
			EdgeDir dir = (o.read[idx].isLeftExtreme()) ? ED_ANTISENSE : ED_SENSE;
			pEdges[idx] = new StringEdge(pVerts[idx], pVerts[1 - idx], dir, comp, overhangs[1 - idx]);
		}
		pEdges[0]->setTwin(pEdges[1]);
		pEdges[1]->setTwin(pEdges[0]);

		pGraph->addEdge(pEdges[0]);
		pGraph->addEdge(pEdges[1]);
	}
	pGraph->writeDot("stringgraph.dot");
	overlapReader.close();
	return pGraph;
}

//
Sequence flip(const Sequence& s)
{
	return reverseComplement(s);
}

//
std::string getOverhangString(const SeqCoord& sc, const std::string& seq)
{
	size_t left, right;

	if(sc.isLeftExtreme())
	{
		 // the end coordinate includes the last base, so the overhang starts at the next pos
		left = sc.interval.end + 1;
		right = seq.length();
	}
	else
	{
		left = 0;
		right = sc.interval.start;
	}
	
	assert(left <= right);
	return seq.substr(left, right - left);
}
