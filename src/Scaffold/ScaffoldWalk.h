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

        // Find the index of the first instance of pvertex in the walk
        // if pvertex is not found, returns -1
        int findVertex(ScaffoldVertex* pVertex) const;

        void print() const;

    private:
        ScaffoldVertex* m_pStartVertex;
        ScaffoldEdgePtrVector m_edges;

};
typedef std::vector<ScaffoldWalk> ScaffoldWalkVector;

#endif
