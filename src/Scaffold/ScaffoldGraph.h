//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGraph - A graph representing long-distance
// relationships between contigs. This shares
// some functionality with Bigraph/StringGraph
// but that implementation is too tuned for low
// memory usage to be generalized.
//
#ifndef SCAFFOLDGRAPH_H
#define SCAFFOLDGRAPH_H

#include "ScaffoldVertex.h"
#include "HashMap.h"

typedef SparseHashMap<VertexID, ScaffoldVertex*> ScaffoldVertexMap;

class ScaffoldGraph
{
    public:
        ScaffoldGraph();
        ~ScaffoldGraph();

        void loadVertices(const std::string& filename);
        void addVertex(ScaffoldVertex* pVertex);
        void addEdge(ScaffoldVertex* pVertex, ScaffoldEdge* pEdge);

        ScaffoldVertex* getVertex(VertexID id) const;

        // Visit each vertex in the graph and call the visit functor object
        template<typename VF>
        bool visit(VF& vf)
        {
            bool modified = false;
            vf.previsit(this);
            ScaffoldVertexMap::const_iterator iter = m_vertices.begin(); 
            for(; iter != m_vertices.end(); ++iter)
            {
                modified = vf.visit(this, iter->second) || modified;
            }
            vf.postvisit(this);
            return modified;
        }

        void writeDot(const std::string& outFile) const;

    private:

        ScaffoldVertexMap m_vertices;

};

#endif
