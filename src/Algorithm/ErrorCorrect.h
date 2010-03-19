//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrect - Methods to identify and correct read errors
//
#ifndef ERRORCORRECT_H
#define ERRORCORRECT_H

#include "Util.h"
#include "Bigraph.h"

namespace ErrorCorrect
{
	// Perform error correction on the given vertex
	std::string correctVertex(Vertex* pVertex, size_t simpleCutoff);

	// simpleCorrect calls a straightforward consensus sequence
	// based on all the reads overlapping pVertex
	std::string simpleCorrect(Vertex* pVertex);

	// trieCorrect builds tries from the overlapping reads
	// to attempt to account for overcollapsed repeats 
	std::string trieCorrect(Vertex* pVertex);
}

#endif
