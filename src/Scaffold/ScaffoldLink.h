//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldLink - Data structure representing
// an inferred relationship between two scaffold 
// components
//
#ifndef SCAFFOLDLINK_H
#define SCAFFOLDLINK_H

#include "GraphCommon.h"
#include "Edge.h"

enum ScaffoldLinkType
{
    SLT_DISTANCEEST,
    SLT_REFERENCE,
    SLT_INFERRED
};

class ScaffoldLink
{
    public:
        ScaffoldLink() {}
        ScaffoldLink(EdgeDir dir, EdgeComp comp, int dist, 
                     double sd, int np, ScaffoldLinkType slt);

        EdgeDir getDir() const;
        EdgeComp getComp() const;

        int distance;
        double stdDev;
        int numPairs;
        ScaffoldLinkType type;
        EdgeData edgeData; // dir/comp member
};

#endif
