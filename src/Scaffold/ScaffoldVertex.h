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
        size_t getNumEdges() const;
        size_t getSeqLen() const;
        double getAStatistic() const;
        bool isRepeat() const;
        VertexID getID() const;
        GraphColor getColor() const;
        std::string getColorString() const;
        ScaffoldVertexClassification getClassification() const;

        //
        void setAStatistic(double v);
        void setClassification(ScaffoldVertexClassification classification);
        void setColor(GraphColor c);

        //
        void deleteEdge(ScaffoldEdge* pEdge);
        void deleteEdgeAndTwin(ScaffoldEdge* pEdge);
        void deleteEdges();
        void deleteEdgesAndTwins(EdgeDir dir);
        void deleteEdgesAndTwins();

        ScaffoldEdge* findEdgeTo(VertexID id, ScaffoldLinkType type) const;
        ScaffoldEdge* findEdgeTo(VertexID id, EdgeDir dir, EdgeComp comp) const;

        ScaffoldEdgePtrVector getEdges();
        ScaffoldEdgePtrVector getEdges(EdgeDir dir);

        // I/O
        void writeDot(std::ostream* pWriter) const;
        void writeEdgesDot(std::ostream* pWriter) const;

    private:
        VertexID m_id;
        size_t m_seqLen;
        double m_AStatistic;
        ScaffoldEdgePtrVector m_edges;
        ScaffoldVertexClassification m_classification;
        GraphColor m_color;
};

#endif
