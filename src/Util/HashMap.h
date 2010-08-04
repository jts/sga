//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HashMap compability defines. Similar to Shaun Jackman's
// compatability code in abyss (see Common/HashMap)
//
//
#ifndef HASHMAP_H
#define HASHMAP_H

#include "config.h"

#if HAVE_UNORDERED_MAP
# include <unordered_map>
# define HashMap std::unordered_map
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#define HashMap std::tr1::unordered_map
#elif HAVE_EXT_HASH_MAP
# undef __DEPRECATED
#include <ext/hash_map>
#define HashMap __gnu_cxx::hash_map
#else
# error No hash map implementation found
#endif

// Define a string hasher if the unordered map was not found
#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP
#elif HAVE_EXT_HASH_MAP
namespace __gnu_cxx                                                                              
{                                                                                             
  template<> struct hash< std::string >                                                       
  {                                                                                           
    size_t operator()( const std::string& x ) const                                           
    {                                                                                         
      return hash< const char* >()( x.c_str() );                                              
    }                                                                                         
  };                                                                                          
}
#endif

// Ensure the sparse hash is available
#if HAVE_GOOGLE_SPARSE_HASH_MAP
#include <google/sparse_hash_map>
#define SparseHashMap google::sparse_hash_map
#else
#error The google sparse hash is required
#endif 

#endif
