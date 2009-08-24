//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Base bidirectional edge class 
//
#include "Edge.h"
#include "Vertex.h"

//
Vertex* Edge::getStart() const 
{ 
	return m_pStart; 
}

//
Vertex* Edge::getEnd() const 
{ 
	return m_pEnd; 
}

//
VertexID Edge::getStartID() const
{ 
	return m_pStart->getID(); 
}

//
VertexID Edge::getEndID() const 
{ 
	return m_pEnd->getID(); 
}

//
EdgeDir Edge::getDir() const 
{ 
	return m_dir; 
}

//
EdgeComp Edge::getComp() const
{ 
	return m_comp; 
}

//
bool Edge::isSelf() const 
{ 
	return m_pStart == m_pEnd; 
}

// 
EdgeDesc Edge::getTwinDesc() const
{
	return EdgeDesc(getStartID(), getTwinDir(), m_comp);
}

// Merge pEdge into this edge by updating the endpoints
void Edge::merge(const Edge* pEdge)
{
	if(pEdge->getComp() == EC_REVERSE)
		flipComp();
	m_pEnd = pEdge->getEnd();
}

// Output
std::ostream& operator<<(std::ostream& out, const Edge& obj)
{
	out << obj.m_pStart->getID() << "," << obj.m_pEnd->getID() << "," << obj.m_dir << "," << obj.m_comp;
	return out;
}

