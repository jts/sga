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
GraphColor Edge::getColor() const
{
	return m_color;
}

// 
EdgeDesc Edge::getTwinDesc() const
{
	return EdgeDesc(getStartID(), getTwinDir(), m_comp);
}

// Join the edge pEdge into this edge, adding to the start
void Edge::join(const Edge* pEdge)
{
	if(pEdge->getComp() == EC_REVERSE)
		flip();

	m_pStart = pEdge->getStart();

	// Now, update the twin of this edge to extend to the twin of pEdge
	m_pTwin->extend(pEdge->getTwin());
}

// Extend this edge by adding pEdge to the end
void Edge::extend(const Edge* pEdge)
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

