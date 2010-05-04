//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// EdgeDesc - A unique description of an edge 
//
#ifndef EDGEDESC_H
#define EDGEDESC_H

#include "GraphCommon.h"

struct EdgeDesc
{
	EdgeDesc(VertexID i, EdgeDir d, EdgeComp c) : id(i), dir(d), comp(c) {}
	VertexID id;
	EdgeDir dir;
	EdgeComp comp;

	// Operators
	bool operator<(const EdgeDesc& obj) const;
	bool operator==(const EdgeDesc& obj) const;
	friend std::ostream& operator<<(std::ostream& out, const EdgeDesc& ed);
};

#endif
