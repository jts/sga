#ifndef ISEQVERTEX_H
#define ISEQVERTEX_H

// Includes
#include <stdio.h>
#include <set>
#include <vector>
#include <ostream>
#include <iostream>
#include <iterator>
#include "Common.h"
#include "Vertex.h"
#include "Util.h"

// Namespace
using namespace std;

// Typedefs
typedef set<Edge> EdgeSet;
typedef vector<Edge> EdgeVec;
typedef EdgeSet::iterator EdgeSetIter;
typedef EdgeVec::iterator EdgeVecIter;

class SeqVertex : public Vertex
{
	public:
		SeqVertex(VertexID id, int kmer, double coverage, string seq) : Vertex(id), m_kmer(kmer), 
												m_overlap( kmer - 1), 
												m_coverage(coverage),
												m_sequence(seq) {}
												
		virtual ~SeqVertex() {}; 

		// Get sequence
		Sequence getSeq() const { return m_sequence; }
		// The cost of traversing this node
		virtual int cost() const { return m_sequence.length() - m_overlap; }
		virtual void merge(const Vertex* pV2, const Edge& e);
	
	private:
		int m_kmer;
		int m_overlap;
		double m_coverage;
		Sequence m_sequence;
};

#endif
