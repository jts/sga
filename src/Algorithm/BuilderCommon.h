///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BuilderCommon -- Common functions and data
// structures for the abstract graph builders
//
#ifndef BUILDER_COMMON_H
#define BUILDER_COMMON_H

#include "Alphabet.h"
#include "GraphCommon.h"
#include "SGUtil.h"

namespace BuilderCommon
{

// Count the number of extensions above the given threshold
size_t countValidExtensions(const AlphaCount64& ac, size_t threshold);

// Make a de Bruijn graph string 
std::string makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction);

// Add a de Bruijn graph edge to the given graph betwee pX and pY.
void addSameStrandDeBruijnEdges(StringGraph* pGraph, const Vertex* pX, const Vertex* pY, EdgeDir direction);

};

#endif
