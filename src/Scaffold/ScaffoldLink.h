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
    SLT_INFERRED,
    SLT_NOTYPE
};

class ScaffoldLink
{
    public:
        ScaffoldLink() {}
        ScaffoldLink(const std::string& id, EdgeDir dir, EdgeComp comp, 
                     int dist, double sd, int np, ScaffoldLinkType slt);

        EdgeDir getDir() const;
        EdgeComp getComp() const;
        char getTypeCode() const;
        static ScaffoldLinkType getType(char tc);

        // IO
        void parse(const std::string& text);
        friend std::ostream& operator<<(std::ostream& out, const ScaffoldLink& link);

        // data 
        std::string endpointID;
        int distance;
        double stdDev;
        int numPairs;
        ScaffoldLinkType type;
        EdgeData edgeData; // dir/comp member
};

#endif
