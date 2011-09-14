///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenProcess - Generate candidate haplotypes
// from an assembly graph for a stream of sites
//
#include "HapgenProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"

//
//
//
HapgenProcess::HapgenProcess(const HapgenParameters& params) : m_parameters(params)
{
}

//
HapgenProcess::~HapgenProcess()
{
}

void HapgenProcess::processSite(const std::string& refName, size_t start, size_t end)
{
    std::cout << "Processing " << refName << " [" << start << " " << end << "]\n";
}

