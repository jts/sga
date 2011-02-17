//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#include "KmerDistribution.h"
#include <assert.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <limits>

KmerDistribution::KmerDistribution()
{

}

void KmerDistribution::add(int kcount)
{
    m_data[kcount]++;
}

double KmerDistribution::getCumulativeProportionLEQ(int n) const
{
    std::vector<int> countVector = toCountVector();
    int64_t sum = 0;
    for(size_t i = 0; i < countVector.size(); ++i)
    {
        sum += countVector[i];
    }


    std::vector<double> cumulativeVector(countVector.size());
    int64_t runningSum = 0;
    for(size_t i = 0; i < cumulativeVector.size(); ++i)
    {
        runningSum += countVector[i];
        cumulativeVector[i] = (double)runningSum / sum;
        //std::cout << "CUMUL " << i << " " << cumulativeVector[i] << "\n";
    }

    assert(n < (int)cumulativeVector.size());
    return cumulativeVector[n];
}


//
int KmerDistribution::findFirstLocalMinimum() const
{
    std::vector<int> countVector = toCountVector();
    if(countVector.empty())
        return -1;

    std::cout << "CV: " << countVector.size() << "\n";
    int prevCount = countVector[1];
    double threshold = 0.75;
    for(size_t i = 2; i < countVector.size(); ++i)
    {
        int currCount = countVector[i];
        double ratio = (double)currCount / prevCount;
        std::cout << i << " " << currCount << " " << ratio << "\n";
        if(ratio > threshold)
            return i;
        prevCount = currCount;
    }
    return -1;
}

int KmerDistribution::getCensoredMode(size_t n) const
{
    std::vector<int> countVector = toCountVector();
    if(countVector.size() < n)
        return -1;
    int modeIdx = -1;
    int modeCount = -1;

    for(; n < countVector.size(); ++n)
    {
        if(countVector[n] > modeCount)
        {
            modeCount = countVector[n];
            modeIdx = n;
        }
    }
    return modeIdx;
}

// Find the boundary of the kmers that are likely erroneous
// We do this by finding the value between 1 and the trusted mode
// that contributes the fewest 
int KmerDistribution::findErrorBoundary() const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cout << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector();
    if(countVector.empty())
        return -1;

    int runningSum = 0;
    double minContrib = std::numeric_limits<double>::max();
    int idx = -1;
    for(int i = 1; i < mode; ++i)
    {
        runningSum += countVector[i];
        double v = (double)countVector[i] / runningSum;
        if(v < minContrib)
        {
            minContrib = v;
            idx = i;
        }
    }
    return idx;
}

// Find the boundary of the kmers that are likely erroneous
// We do this by finding the value between 1 and the trusted mode
// that contributes the fewest 
int KmerDistribution::findErrorBoundaryByRatio(double ratio) const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cout << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector();
    if(countVector.empty())
        return -1;

    for(int i = 1; i < mode - 1; ++i)
    {
        int currCount = countVector[i];
        int nextCount  = countVector[i+1];
        double cr = (double)currCount / nextCount;
        if(cr < ratio)
            return i;
    }
    return -1;
}

// 
std::vector<int> KmerDistribution::toCountVector() const
{
    std::vector<int> out;
    if(m_data.empty())
        return out;

    int min = 0;
    int max = 200;//m_data.rbegin()->first;

    for(int i = min; i <= max; ++i)
    {
        std::map<int,int>::const_iterator iter = m_data.find(i);
        int v = (iter != m_data.end()) ? iter->second : 0;
        out.push_back(v);
    }
    return out;
}

void KmerDistribution::print(int max) const
{
    printf("Kmer coverage histogram\n");
    printf("cov\tcount\n");

    int maxCount = 0;
    std::map<int,int>::const_iterator iter = m_data.begin();
    for(; iter != m_data.end(); ++iter)
    {
        if(iter->first <= max)
            printf("%d\t%d\n", iter->first, iter->second);
        else
            maxCount += iter->second;
    }
    printf(">%d\t%d\n", max, maxCount);

}
