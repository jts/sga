//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#ifndef KMERDISTRIBUTION_H
#define KMERDISTRIBUTION_H

#include <vector>
#include <map>

class KmerDistribution
{
    public:
    
        KmerDistribution();
        ~KmerDistribution() {}

        void add(int count);
        void print(int max) const; 

    private:

        // int -> int map of the number of times a kmer with multiplicty N has been seen
        std::map<int,int> m_data;
};

#endif
