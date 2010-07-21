//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#ifndef VERTEX_H
#define VERTEX_H

// Includes
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"
#include "TransitiveGroupCollection.h"
#include "QualityVector.h"
#include "EncodedString.h"
#include "SimpleAllocator.h"

// Forward declare
class Edge;

// Default edge sorting function, by ID
struct EdgeIDComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};

// Edge sorting function, by length
struct EdgeLenComp
{
    bool operator()(const Edge* pA, const Edge* pB);
};


// Typedefs
typedef std::map<EdgeDesc, Edge*> EdgePtrMap;
typedef std::vector<Edge*> EdgePtrVec;
typedef std::set<EdgeDesc> EdgeDescSet;
typedef std::list<Edge*> EdgePtrList;
typedef EdgePtrMap::iterator EdgePtrMapIter;
typedef EdgePtrMap::const_iterator EdgePtrMapConstIter;
typedef EdgePtrVec::iterator EdgePtrVecIter;
typedef EdgePtrVec::const_iterator EdgePtrVecConstIter;
typedef EdgePtrList::iterator EdgePtrListIter;
typedef EdgePtrList::const_iterator EdgePtrListConstIter;

class Vertex
{
    public:
    
        Vertex(VertexID id, const std::string& s) : m_id(id), 
                                                    m_seq(s), 
                                                    m_color(GC_WHITE),
                                                    m_isContained(false) {}
        ~Vertex();

        // High-level modification functions
        
        // Merge another vertex into this vertex, as specified by pEdge
        void merge(Edge* pEdge);

        // sort the edges by the ID of the vertex they point to
        void sortAdjListByID();

        // sort the edges by the length of the label of the edge
        void sortAdjListByLen();

        // Ensure that all the edges are unique
        void makeUnique(); 

        // Compute the transitive groups for this vertex
        // A transitive group is a set of edges s.t. one
        // edge is irreducible and the other edges in the set
        // are transitive w.r.t. the irreducible edge.
        TransitiveGroupCollection computeTransitiveGroups(EdgeDir dir);

        // Get a multioverlap object representing the overlaps for this vertex
        MultiOverlap getMultiOverlap() const;

        // Construct a trie from the edges, one for each each direction
        void fillTries(double p_error, SeqTrie* pSenseTrie, SeqTrie* pAntisenseTrie) const;

        // Return the inferred quality value for each base in the sequence
        // using the overlap information
        QualityVector getInferredQuality() const;

        //
        std::string getInferredConsensus() const;

        // Return the prior probabiltiy based on quality scores
        QualityVector getPriorQuality() const;

        // Edge list operations
        void addEdge(Edge* ep);
        void removeEdge(Edge* pEdge);
        void removeEdge(const EdgeDesc& ed);
        void deleteEdge(Edge* pEdge);
        void deleteEdges();
        void sweepEdges(GraphColor c);
        bool hasEdge(Edge* pEdge) const;
        bool hasEdge(const EdgeDesc& ed) const;
        bool hasEdgeTo(const Vertex* pY) const;

        Edge* getEdge(const EdgeDesc& ed);
        EdgePtrVec findEdgesTo(VertexID id);
        EdgePtrVec getEdges(EdgeDir dir) const;
        EdgePtrVec getEdges() const;
        EdgePtrVecIter findEdge(const EdgeDesc& ed);
        EdgePtrVecConstIter findEdge(const EdgeDesc& ed) const;

        size_t countEdges() const;
        size_t countEdges(EdgeDir dir);

        // Ensure the vertex data is sane
        void validate() const;
        
        // setters
        void setID(VertexID id) { m_id = id; }
        void setEdgeColors(GraphColor c);
        void setSeq(const std::string& s) { m_seq = s; }
        void setColor(GraphColor c) { m_color = c; }
        void setContained(bool c) { m_isContained = c; }

        // getters
        VertexID getID() const { return m_id; }
        GraphColor getColor() const { return m_color; }
        const DNAEncodedString& getSeq() const { return m_seq; }
        std::string getStr() const { return m_seq.toString(); }
        size_t getSeqLen() const { return m_seq.length(); }
        size_t getMemSize() const;
        bool isContained() const { return m_isContained; }

        // Memory management
        void* operator new(size_t /*size*/)
        {
            return SimpleAllocator<Vertex>::Instance()->alloc();
        }

        void operator delete(void* target, size_t /*size*/)
        {
            SimpleAllocator<Vertex>::Instance()->dealloc(target);
        }

        // Output edges in graphviz format
        void writeEdges(std::ostream& out, int dotFlags) const;

    private:

        // Ensure all the edges in DIR are unique
        void makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec);

        VertexID m_id;
        EdgePtrVec m_edges;
        DNAEncodedString m_seq;
        GraphColor m_color;
        bool m_isContained;
};

#endif
