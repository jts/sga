//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldEdge - An edge in a scaffold graph. 
//
#include "ScaffoldEdge.h"

ScaffoldEdge::ScaffoldEdge(ScaffoldVertex* pEnd, Edge* pTwin, 
                           EdgeDir dir, EdgeComp comp) : m_pEnd(pEnd), m_pTwin(pTwin)
{
    m_edgeData.setDir(dir);
    m_edgeData.setComp(comp);
}

