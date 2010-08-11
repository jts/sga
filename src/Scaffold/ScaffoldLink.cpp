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
#include "ScaffoldLink.h"

ScaffoldLink::ScaffoldLink(EdgeDir dir, EdgeComp comp, 
                           int dist, double sd, 
                           int np, ScaffoldLinkType slt) : distance(dist), stdDev(sd), 
                                                           numPairs(np), type(slt)
{
    edgeData.setDir(dir);
    edgeData.setComp(comp);
}

//
EdgeDir ScaffoldLink::getDir() const
{
    return edgeData.getDir();
}

//
EdgeComp ScaffoldLink::getComp() const
{
    return edgeData.getComp();
}
