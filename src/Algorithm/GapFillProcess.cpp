///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapFillProcess - Fill in gaps in a scaffold
//
#include <algorithm>
#include "GapFillProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "HaplotypeBuilder.h"
#include "MultiAlignment.h"

//
//
//
GapFillProcess::GapFillProcess(const GapFillParameters& params) : m_parameters(params)
{
}

//
GapFillProcess::~GapFillProcess()
{
}

void GapFillProcess::processScaffold(const std::string& scaffold)
{
    std::cout << "Processing scaffold of length " << scaffold.length() << "\n";
}

