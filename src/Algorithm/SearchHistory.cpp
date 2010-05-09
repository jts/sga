//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// SearchHistory - data structures to track
// the history of a FM-index search
//
//-----------------------------------------------
#include "SearchHistory.h"
#include <algorithm>
#include <iterator>

//
// Link
//
SearchHistoryLink::SearchHistoryLink() : pNode(NULL)
{

}

//
SearchHistoryLink::SearchHistoryLink(SearchHistoryNode* ptr) : pNode(ptr)
{
    if(pNode != NULL)
        pNode->increment();
}

// We need to handle the copy constructor and the assignment operator
// for the reference counting to be correct
SearchHistoryLink::SearchHistoryLink(const SearchHistoryLink& link) : pNode(link.pNode)
{
    if(pNode != NULL)
        pNode->increment();
}

//
SearchHistoryLink& SearchHistoryLink::operator=(SearchHistoryLink const& link)
{
    SearchHistoryNode* pOld = pNode;
    pNode = link.pNode;
    if(pNode != NULL)
        pNode->increment();

    if(pOld != NULL)
    {
        pOld->decrement();
        if(pOld->getCount() == 0)
            delete pOld;
    }
    return *this;
}

//
SearchHistoryLink::~SearchHistoryLink()
{
    if(pNode != NULL)
    {
        pNode->decrement();
        if(pNode->getCount() == 0)
            delete pNode;
    }
}


//
// Node
//


// adding a child of the node automatically increments the refCount of this node
// through the child's Link to this node. Once all the children of this node
// have been removed it will automatically be deleted
SearchHistoryLink SearchHistoryNode::createChild(int var_pos, char var_base)
{
    return SearchHistoryLink(new SearchHistoryNode(this, var_pos, var_base));
}

// The root has NULL as a parent 
SearchHistoryLink SearchHistoryNode::createRoot()
{
    return SearchHistoryLink(new SearchHistoryNode(NULL, -1, ROOT_CHAR));
}

// Return the search history up to the root node
SearchHistoryVector SearchHistoryNode::getHistoryVector()
{
    SearchHistoryVector out;
    SearchHistoryLink pCurr(this);
    bool done = false;
    while(!done)
    {
        if(pCurr->m_variant.base == ROOT_CHAR)
        {
            // root found, terminate the search
            done = true;
        }
        else
        {
            out.add(pCurr->m_variant);
            pCurr = pCurr->m_parentLink;
        }
    }
    return out;
}

//
// SearchHistoryVector
//

// Normalize the search history so that the positions are in ascending (left to right)
// order and the bases are wrt. the original strand of the forward sequence
// This is required to compare two SearchHistories that represent alignments to different 
// strands
void SearchHistoryVector::normalize(bool doComplement)
{
    std::sort(m_history.begin(), m_history.end(), SearchHistoryItem::sortPos);

    if(doComplement)
    {
        for(size_t i = 0; i < m_history.size(); ++i)
        {
            m_history[i].base = complement(m_history[i].base);    
        }
    }
}

//
void SearchHistoryVector::add(int pos, char base)
{
    SearchHistoryItem item(pos, base);
    m_history.push_back(item);
}

// 
void SearchHistoryVector::add(SearchHistoryItem& item)
{
    m_history.push_back(item);
}

//
int SearchHistoryVector::countDifferences(const SearchHistoryVector& a, const SearchHistoryVector& b, int maxPos)
{
    size_t na = a.m_history.size();
    size_t nb = b.m_history.size();

    size_t i = 0;
    size_t j = 0;
    int count = 0;

    while(i < na && j < nb)
    {
        const SearchHistoryItem& itemA = a.m_history[i];
        const SearchHistoryItem& itemB = b.m_history[j];
        
        if(itemA.pos > maxPos || itemB.pos > maxPos)
            break;

        if(itemA.pos == itemB.pos)
        {
            if(itemA.base != itemB.base)
                ++count;
            ++i;
            ++j;
        }
        else if(itemA.pos < itemB.pos)
        {
            ++count;
            ++i;
        }
        else if(itemB.pos < itemA.pos)
        {
            ++count;
            ++j;
        }
        else
        {
            assert(false);
        }
    }
    
    // Count the remaining elements in A and B that are less than maxpos
    while(i < na && a.m_history[i].pos <= maxPos)
    {
        ++count;
        ++i;
    }

    // Count the remaining elements in A and B that are less than maxpos
    while(j < nb && b.m_history[j].pos <= maxPos)
    {
        ++count;
        ++j;
    }

    return count;
}

//
std::ostream& operator<<(std::ostream& out, const SearchHistoryVector& hist)
{
    std::copy(hist.m_history.begin(), hist.m_history.end(), 
              std::ostream_iterator<SearchHistoryItem>(out, "\t"));
    return out;
}

