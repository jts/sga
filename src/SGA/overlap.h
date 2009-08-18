//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// overlap - Overlap reads using a bwt
//
#ifndef OVERLAP_H
#define OVERLAP_H

void printOverlapUsage();
int overlapMain(int argc, char** argv);
void computeOverlaps(std::string indexFile, std::string overlapFile);

#endif
