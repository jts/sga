//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiAlignment - Class constructing a multiple
// alignment between a root sequence and a set of sequences
//
#include <iostream>
#include <algorithm>
#include "MultiAlignment.h"
#include "Util.h"

struct CigarIter
{
    const MAlignData* pData;
    int symIdx; // the index of the cigar string we are parsing
    int baseIdx; // the index of the current base 

    char updateAndEmit(char updateMode)
    {
        char cigSym = getCigarSymbol();
        char outBase;

        if(updateMode == 'I')
        {
            if(cigSym == 'I')
            {
                symIdx += 1;
                outBase = '-';
            }
            else if(cigSym == 'D')
            {
                assert(false);
            }
            else
            {
                symIdx += 1;
                outBase = getOutputBase(cigSym);
                if(outBase != '.' && outBase != '-')
                    baseIdx += 1;
            }
        }
        else if(updateMode == 'D')
        {
            if(cigSym == 'D')
            {
                symIdx += 1;
                outBase = getOutputBase(cigSym);
                if(outBase != '.')
                    baseIdx += 1;
            }
            else if(cigSym == 'I')
            {
                assert(false);
            }
            else
            {
                outBase = '-';
                //symIdx += 1;
            }
        }
        else
        {
            symIdx += 1;
            outBase = getOutputBase(cigSym);
            if(outBase != '.')
                baseIdx += 1;
        }

        return outBase;
    }

    char getCigarSymbol()
    {
        if(symIdx >= (int)pData->expandedCigar.size())
            return 'S';
        return pData->expandedCigar[symIdx];
    }

    char getOutputBase(char cigSym)
    {
        if(cigSym == 'S')
            return '.';
        else
            return pData->str[baseIdx];
    }

    static bool sortPosition(const CigarIter& a, const CigarIter& b)
    {
        return a.pData->position < b.pData->position;
    }
};
typedef std::vector<CigarIter> CigarIterVector;

MultiAlignment::MultiAlignment(std::string rootStr, const MAlignDataVector& inData)
{
    // Build a padded multiple alignment from the pairwise alignments to the root
    CigarIterVector iterVec;

    // Create entry for the root
    MAlignData rootData = {rootStr, std::string(rootStr.size(), 'M'), 0};
    CigarIter tmp = {&rootData, 0, 0};
    iterVec.push_back(tmp);
    printf("%zu\t%s\n", 0, tmp.pData->expandedCigar.c_str());
    for(size_t i = 0; i < inData.size(); ++i)
    {
        CigarIter iter = {&inData[i], 0, 0};
        iterVec.push_back(iter);
        printf("%zu\t%s\n", i+1, iter.pData->expandedCigar.c_str());
    }

    std::stable_sort(iterVec.begin(), iterVec.end(), CigarIter::sortPosition);

    bool done = false;
    StringVector outStrings(iterVec.size());
    while(!done)
    {
        // Check if any strings have a deletion or insertion at this base
        bool hasDel = false;
        bool hasInsert = false;
        bool hasMatch = false;

        for(size_t i = 0; i < iterVec.size(); ++i)
        {
            char sym = iterVec[i].getCigarSymbol();
            if(sym == 'D')
                hasDel = true;
            else if(sym == 'I')
                hasInsert = true;
            else if(sym == 'M')
                hasMatch = true;
        }

        done = !hasDel && !hasInsert && !hasMatch;
        if(done)
            break;

        char updateMode = hasDel ? 'D' : (hasInsert ? 'I' : 'M');
        for(size_t i = 0; i < iterVec.size(); ++i)
        {
            char outSym = iterVec[i].updateAndEmit(updateMode);
            outStrings[i].append(1,outSym);
        }
    }

    size_t len = outStrings[0].size();
    for(size_t l = 0; l < len; l += 80)
    {
        for(size_t i = 0; i < outStrings.size(); ++i)
        {
            size_t stop = outStrings[i].size() - l < 80 ? outStrings[i].size() - l : 80;        
            printf("%zu\t%s\n", i, outStrings[i].substr(l, stop).c_str());
        }

        std::cout << "\n";
    }
}

