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
// Called from Edge::flip
void SEData::flip()
{
	seq = reverseComplement(seq);
}

// Get the edge's label (for writeDot)
std::string SEData::getLabel() const
{
	return seq;
}

// Implementation of edge merge
void SVData::merge(const SVData& other, const StringEdge& edge)
{
	// Merge the sequence
	if(edge.getDir() == ED_SENSE)
		seq.append(edge.getData().seq);
	else
		seq.insert(0, edge.getData().seq);

	// This is only valid for perfect matches
	std::cout << "this:  " << seq << "\n"; 
	std::cout << "other: " << other.seq << "\n";
	assert(seq.find(other.seq) != std::string::npos);
	readCount += other.readCount;
}

// Visitor which outputs the graph in fasta format
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, StringVertex* pVertex)
{
	const SVData& data = pVertex->getData();
	m_fileHandle << ">" << pVertex->getID() << " " << data.readCount << "\n";
	m_fileHandle << data.seq << "\n";
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
		pGraph->addVertex(new StringVertex(si.id, SVData(si.id, si.seq)));

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
			pVerts[idx] = pGraph->getVertex(o.read[idx].id);
			assert(pVerts[idx]);
			sizes[idx] = pVerts[idx]->getData().seq.length();
			
			// Ensure that the overlaps are canonical
			assert(o.read[idx].interval.start < o.read[idx].interval.end);
			
			// Ensure the reads are not identical
			assert(!o.read[idx].isContainment(sizes[idx]));

			overhangs[idx] = getOverhangString(o.read[idx], pVerts[idx]->getData().seq);
		}

		// Ensure that each overlap is extreme
		if(!o.read[0].isExtreme(sizes[0]) || !o.read[1].isExtreme(sizes[1]))
			continue;

		EdgeComp comp = (o.read[0].isReverse() || o.read[1].isReverse()) ? EC_REVERSE : EC_SAME;

		// Add edges
		for(size_t idx = 0; idx < 2; ++idx)
		{
			EdgeDir dir = (o.read[idx].isLeftExtreme()) ? ED_ANTISENSE : ED_SENSE;
			StringEdge e(o.read[idx].id, o.read[1 - idx].id, dir, comp, SEData(overhangs[1 - idx]));
			pGraph->addEdge(e);
		}
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
