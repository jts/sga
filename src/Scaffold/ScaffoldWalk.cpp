//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldWalk - a walk through a scaffold graph
//
#include "ScaffoldWalk.h"

ScaffoldWalk::ScaffoldWalk(const ScaffoldEdgePtrVector& edges)
{
    m_edges = edges;
}

void ScaffoldWalk::print() const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter)
    {
        std::cout << *(*iter) << " ";
    }
    std::cout << "\n";
}
