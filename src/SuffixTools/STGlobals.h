//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// STGlobals.h - Common definitions and declarations for suffix tree,
// suffix array and BWT data types
//
#ifndef STGLOBALS_H
#define STGLOBALS_H
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
#include <set>
#include <assert.h>

//
// Constants
//
const uint8_t ALPHABET_SIZE = 5;
const char ALPHABET[ALPHABET_SIZE] = {'A', 'C', 'G', 'T', '$'};

//
// Basic typedefs
//
typedef std::string BWStr;
typedef uint32_t Label;
typedef uint8_t AIdx;
typedef std::vector<int> IntVector;
typedef std::vector<std::string> StringVector;

#endif
