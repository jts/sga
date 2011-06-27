///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#include "GraphCompare.h"
#include "BWTAlgorithms.h"

GraphCompare::GraphCompare(const BWT* pBaseBWT, 
                           const BWT* pBaseRBWT,
                           const BWT* pVariantBWT, 
                           const BWT* pVariantRBWT, 
                           int kmer) : m_pBaseBWT(pBaseBWT),  
                                       m_pBaseRevBWT(pBaseRBWT),
                                       m_pVariantBWT(pVariantBWT),  
                                       m_pVariantRevBWT(pVariantRBWT),
                                       m_kmer(kmer)
{


}

//
GraphCompare::~GraphCompare()
{

}

// Run the actual comparison
void GraphCompare::run()
{
    std::cout << "Running graph comparison\n";
}   
