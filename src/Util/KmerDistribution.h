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
#include <cstddef>

class KmerDistribution
{
    public:
    
        KmerDistribution();
        ~KmerDistribution() {}
        
        // Returns the proportion of the distribution less than or equal to n
        double getCumulativeProportionLEQ(int n) const;

        // Returns the predicted boundary of the left tail which holds erroneous kmers
        int findErrorBoundary() const;
        int findErrorBoundaryByRatio(double ratio) const;

        // Get the mean value of the distribution if the first n values are ignored
        int getCensoredMode(size_t n) const;

        //
        int findFirstLocalMinimum() const;
        void add(int count);
        void print(int max) const; 

    private:
        std::vector<int> toCountVector() const;

        // int -> int map of the number of times a kmer with multiplicty N has been seen
        std::map<int,int> m_data;
};

#endif
