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
        SeqVertex(VertexID id, string seq) : Vertex(id), m_sequence(seq) {}
        virtual ~SeqVertex() {}; 

        // Get sequence
        Sequence getSeq() const { return m_sequence; }
        virtual void merge(const Vertex* pV2, const Edge& e);
    
    private:
        Sequence m_sequence;
};

#endif
