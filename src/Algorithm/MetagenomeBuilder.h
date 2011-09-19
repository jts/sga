///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MetagenomeBuilder -- Implementation of the metagenome
// assembly process. This class takes in a kmer and 
// assembles a contig locally around that sequence
//
#ifndef METAGENOME_BUILDER_H
#define METAGENOME_BUILDER_H

#include "BWT.h"
#include "SGUtil.h"
#include "VariationBubbleBuilder.h"

class MetagenomeBuilder
{
    public:

        MetagenomeBuilder();
        ~MetagenomeBuilder();
        
        // Set parameters
        void setSource(const std::string& seq, int coverage);
        void setIndex(const BWT* pBWT, const BWT* pRevBWT);
        void setKmerParameters(size_t k, size_t threshold);

        // run the assembly
        void run();
        
        // Get the contigs from the graph
        void getContigs(StringVector& contigs);

    private:
        
        // Functions
        // Add a vertex to the graph
        void addVertex(Vertex* pVertex, int coverage);

        // Data

        StringGraph* m_pGraph;
        StrIntMap m_vertexCoverageMap;
        BuilderExtensionQueue m_queue;

        const BWT* m_pBWT;
        const BWT* m_pRevBWT;
        size_t m_kmer;
        size_t m_kmerThreshold;
};

#endif
