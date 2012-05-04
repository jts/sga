///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VariationBuilderCommon - Data structures
// commons to various algorithms for assembling
// variation haplotypes
//
#ifndef VARIATION_BUILDER_COMMON_H
#define VARIATION_BUILDER_COMMON_H
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include <queue>

// Result of the bubble construction
// Used to track statistics
enum BubbleResultCode
{
    BRC_UNKNOWN,
    BRC_OK,
    BRC_SOURCE_BROKEN,
    BRC_SOURCE_BRANCH,
    BRC_TARGET_BROKEN,
    BRC_TARGET_BRANCH,
    BRC_WALK_FAILED,
    BRC_HB_FAILED,
    BRC_NO_SOLUTION
};

// The actual result structure
struct BubbleResult
{
    std::string targetString;
    std::string sourceString;
    double targetCoverage;
    double sourceCoverage;
    BubbleResultCode returnCode;
};

// A directed node in the de Bruijn graph that
// has not had it's neighbors visited
struct BuilderExtensionNode
{
    BuilderExtensionNode(Vertex* pX, EdgeDir d) : pVertex(pX), direction(d), distance(0) {}
    BuilderExtensionNode(Vertex* pX, EdgeDir dir, int dist) : pVertex(pX), direction(dir), distance(dist) {}

    Vertex* pVertex; // the vertex to extend
    EdgeDir direction; // the direction to extend to
    int distance; // the total number of nodes from the start to this node.
};
typedef std::queue<BuilderExtensionNode> BuilderExtensionQueue;
typedef std::map<std::string, int> StrIntMap;

#endif
