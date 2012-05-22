//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PopulationIndex - A structure mapping a read
// index to a member of a population
//
#ifndef POPULATION_INDEX_H
#define POPULATION_INDEX_H

#include <vector>
#include <string>

struct PopulationMember
{
    size_t start;
    size_t end;
    std::string name;

    //
    static bool compareByStart(const PopulationMember& a,
                               const PopulationMember& b)
    {
        return a.start < b.start;
    }
};

class PopulationIndex
{
    public:
        PopulationIndex(const std::string& filename);

        // Return the name of the individual for the read at index
        std::string getName(size_t index) const;

    private:
        std::vector<PopulationMember> m_population;
};

#endif
