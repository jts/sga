//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSearch - Find walks through a scaffold graph
//
#ifndef SCAFFOLDSEARCH_H
#define SCAFFOLDSEARCH_H

#include "GraphSearchTree.h"
#include "ScaffoldWalk.h"

typedef std::vector<ScaffoldEdgePtrVector> SEPVVector;

namespace ScaffoldSearch
{
    void findVariantWalks(ScaffoldVertex* pX, 
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          ScaffoldWalkVector& outWalks);

    void convertEdgeVectorsToScaffoldWalk(const SEPVVector& edgeWalks, 
                                          ScaffoldWalkVector& outWalks);
    
};

#endif
