//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldEdge - An edge in a scaffold graph. 
//
#ifndef SCAFFOLDEDGE_H
#define SCAFFOLDEDGE_H

#include "Edge.h"

class ScaffoldVertex;

class ScaffoldEdge
{
    public:

        ScaffoldEdge(ScaffoldVertex* pEnd, Edge* pTwin, EdgeDir dir, EdgeComp comp);

    private:
        ScaffoldVertex* m_pEnd;
        Edge* m_pTwin;
        EdgeData m_edgeData; // dir/comp member

};

#endif
