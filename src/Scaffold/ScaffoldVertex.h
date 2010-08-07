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

enum ScaffoldVertexClassification
{
    SVC_UNIQUE,
    SVC_REPEAT,
    SVC_UNKNOWN
};

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
        std::string getColorString() const;
        double getAStatistic() const;
        ScaffoldVertexClassification getClassification() const;

        //
        void setAStatistic(double v);
        void setClassification(ScaffoldVertexClassification classification);

        //
        ScaffoldEdge* findEdgeTo(VertexID id, ScaffoldEdgeType type);
            
        // I/O
        void writeDot(std::ostream* pWriter) const;
        void writeEdgesDot(std::ostream* pWriter) const;

    private:
        VertexID m_id;
        size_t m_seqLen;
        double m_AStatistic;
        ScaffoldEdgePtrVector m_edges;
        ScaffoldVertexClassification m_classification;
};

#endif
