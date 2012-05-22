//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PopulationIndex - A structure mapping a read
// index to a member of a population
//
#include "PopulationIndex.h"
#include "Util.h"
#include <algorithm>

PopulationIndex::PopulationIndex(const std::string& filename)
{
    std::cout << "Loading population index from " << filename << "\n";

    std::istream* reader = createReader(filename);
    std::string line;
    while(getline(*reader, line))
    {
        std::stringstream parser(line);
        PopulationMember member;
        parser >> member.start;
        parser >> member.end;
        parser >> member.name;
        m_population.push_back(member);
    }
    delete reader;
    reader = NULL;

    // Make sure the index is properly formatted
    assert(m_population.size() != 0);
    assert(m_population[0].start == 0);
    for(size_t i = 1; i < m_population.size(); ++i)
    {
        assert(m_population[i].start > m_population[i-1].start);
        assert(m_population[i].start == m_population[i-1].end + 1);
    }
    
    /* Test the index
    size_t start = 0;
    size_t end = m_population.back().end + 2;
    for(size_t i = start; i <= end; ++i)
    {
        std::string name = getName(i);
        printf("%zu = %s\n", i, name.c_str());
    }
    */
}

//
std::string PopulationIndex::getName(size_t index) const
{
    assert(index <= m_population.back().end);
    // Perform a binary search to identify the first element of the population
    // that is strictly greater than the index. This returns an iterator that is 
    // guarenteed to be one past the individual that index belongs to.
    PopulationMember search = { index, index, "" };
    std::vector<PopulationMember>::const_iterator iter = std::upper_bound(m_population.begin(), 
                                                                          m_population.end(), 
                                                                          search, 
                                                                          PopulationMember::compareByStart);
    assert(iter != m_population.begin());
    iter--;
    assert(index >= iter->start && index <= iter->end);
    return iter->name;
}

