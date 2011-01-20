//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGroup - A group of links that
// can be potentially ordered into a scaffold 
//
#ifndef SCAFFOLDGROUP_H
#define SCAFFOLDGROUP_H

#include "ScaffoldVertex.h"
#include "ScaffoldLink.h"

typedef std::list<ScaffoldLink> LinkList;
typedef LinkList::iterator LinkListIterator;

class ScaffoldGroup
{
    public:
        ScaffoldGroup(const ScaffoldVertex* pRootVertex, int maxOverlap);

        void addLink(const ScaffoldLink& link);
        void computeBestOrdering();
        int scoreLinkPlacement(const ScaffoldLink& link, const LinkList& unplacedList);
        int calculateLongestOverlap();

    private:

        const ScaffoldVertex* m_pRootVertex;
        int m_maxOverlap;
        std::vector<ScaffoldLink> m_links;

};

#endif
