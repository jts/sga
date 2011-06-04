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
#include <map>
#include "MultiAlignment.h"

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
                outBase = '-';
            }
            else if(cigSym == 'S')
            {
                outBase = '.';
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
    m_verbose = 1;
    // Build a padded multiple alignment from the pairwise alignments to the root
    CigarIterVector iterVec;

    // Create entry for the root
    MAlignData rootData;
    rootData.str = rootStr;
    rootData.expandedCigar = std::string(rootStr.size(), 'M');
    rootData.position = 0;
    rootData.targetID = -1;
    rootData.targetAlignLength = - 1;
    rootData.targetPosition = -1;

    m_alignData.push_back(rootData);
    m_alignData.insert(m_alignData.end(), inData.begin(), inData.end());

    CigarIter tmp = {&rootData, 0, 0};
    iterVec.push_back(tmp);

    if(m_verbose > 1)
        printf("%d\t%s\n", 0, tmp.pData->expandedCigar.c_str());

    for(size_t i = 0; i < inData.size(); ++i)
    {
        CigarIter iter = {&inData[i], 0, 0};
        iterVec.push_back(iter);
        if(m_verbose > 1)
            printf("%zu\t%s\n", i+1, iter.pData->expandedCigar.c_str());
    }

//    std::stable_sort(iterVec.begin(), iterVec.end(), CigarIter::sortPosition);

    bool done = false;
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
            m_alignData[i].padded.append(1,outSym);
        }
    }

    generateConsensus();
}

// Generate simple consensus string from the multiple alignment
std::string MultiAlignment::generateConsensus()
{
    assert(!m_alignData.empty() && !m_alignData.front().padded.empty());

    std::string consensus;
    std::string paddedConsensus;

    std::map<char, int> baseMap;

    size_t num_rows = m_alignData.size() - 1; // do not include root
    size_t num_cols = m_alignData.front().padded.size();
    for(size_t i = 0; i < num_cols; ++i)
    {
        baseMap.clear();
        for(size_t j = 1; j < num_rows; ++j)
        {
            char b = m_alignData[j].padded[i];
            if(b != '.')
                baseMap[b]++;
        }

        int max = 0;
        char maxBase = '.';
        for(std::map<char,int>::iterator iter = baseMap.begin(); iter != baseMap.end(); ++iter)
        {
            if(iter->second > max)
            {
                max = iter->second;
                maxBase = iter->first;
            }
        }

        paddedConsensus.append(1,maxBase);

        if(maxBase == '.' && consensus.empty())
            continue; // skip no-call at beginning
        else if(maxBase == '.')
            break; // non-called position in middle of read, stop consensus generation
        else if(maxBase != '-') //
            consensus.append(1, maxBase);
    }

    if(m_verbose > 0)
        print(&paddedConsensus);

    return consensus;
}

//
void MultiAlignment::print(const std::string* pConsensus) const
{
    assert(!m_alignData.empty() && !m_alignData.front().padded.empty());
    size_t len = m_alignData[0].padded.size();
    int col_size = 120;
    for(size_t l = 0; l < len; l += col_size)
    {
        if(pConsensus != NULL)
        {
            int diff = pConsensus->size() - l;
            int stop = diff < col_size ? diff : col_size;
            if(stop > 0)
                printf("C\t%s\n", pConsensus->substr(l,stop).c_str());
            else
                printf("C\n");
        }
        
        for(size_t i = 0; i < m_alignData.size(); ++i)
        {
            const MAlignData& mad = m_alignData[i];
            int diff = mad.padded.size() - l;
            int stop = diff < col_size ? diff : col_size;
            printf("%zu\t%s\t%d,%d,%d\n", i, mad.padded.substr(l, stop).c_str(), mad.targetID, mad.targetPosition, mad.targetAlignLength);
        }

        std::cout << "\n";
    }
}
