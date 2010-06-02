//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// TransitiveGroupCollection - A collection
// of TransitiveGroups, representing the partitioning
// of edges for a given vertex
//
#ifndef TRANSITIVEGROUPCOLLECTION_H
#define TRANSITIVEGROUPCOLLECTION_H

#include "GraphCommon.h"
#include "TransitiveGroup.h"

class Vertex;
struct EdgeDesc;

typedef std::vector<TransitiveGroup> TransitiveGroupVector;

class TransitiveGroupCollection
{
    public:
        TransitiveGroupCollection(Vertex* pVertex, EdgeDir dir);

        TransitiveGroup& createGroup(Edge* pIrreducible);
        TransitiveGroup& getGroup(size_t idx);
        size_t findGroup(EdgeDesc ed) const;
        size_t numGroups() const;

        void print() const;

    private:

        Vertex* m_pVertex;
        EdgeDir m_dir;
        TransitiveGroupVector m_groups;
};

#endif
