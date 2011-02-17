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

KmerDistribution::KmerDistribution()
{

}

void KmerDistribution::add(int kcount)
{
    m_data[kcount]++;
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
