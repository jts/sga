//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldWalk - a walk through a scaffold graph
//
#ifndef SCAFFOLDWALK_H
#define SCAFFOLDWALK_H

#include <vector>
#include "ScaffoldEdge.h"

typedef std::vector<ScaffoldEdge*> ScaffoldEdgePtrVector;

// A walk through the scaffold graph
class ScaffoldWalk
{
    public:
        ScaffoldWalk(ScaffoldVertex* pStartVertex);
        void addEdge(ScaffoldEdge* pEdge);

        void print() const;

    private:
        ScaffoldEdgePtrVector m_edges;

};
typedef std::vector<ScaffoldWalk> ScaffoldWalkVector;

#endif
