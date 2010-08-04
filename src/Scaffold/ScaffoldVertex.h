//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVertex - A node in a scaffold graph. 
//
#ifndef SCAFFOLDVERTEX_H
#define SCAFFOLDVERTEX_H

#include "Bigraph.h"
#include "ScaffoldEdge.h"

typedef std::vector<ScaffoldEdge*> ScaffoldEdgePtrVector;
class ScaffoldVertex
{
    public:
        
        //
        ScaffoldVertex(VertexID id, size_t seqLen);
        ~ScaffoldVertex();
        
        //
        void addEdge(ScaffoldEdge* pEdge);
        
        //
        VertexID getID() const;
        size_t getNumEdges() const;
        size_t getSeqLen() const;

        //
        ScaffoldEdge* findEdgeTo(VertexID id, ScaffoldEdgeType type);
            
        // I/O
        void writeEdgesDot(std::ostream* pWriter) const;

    private:
        VertexID m_id;
        size_t m_seqLen;
        ScaffoldEdgePtrVector m_edges;

};

#endif
