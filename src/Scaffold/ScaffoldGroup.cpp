//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGroup - A set of links that
// can be potentially ordered into a scaffold 
//
#include "ScaffoldGroup.h"
#include "Interval.h"

//
ScaffoldGroup::ScaffoldGroup(const ScaffoldVertex* pRootVertex, 
                             int maxOverlap) : m_pRootVertex(pRootVertex),
                                               m_maxOverlap(maxOverlap)
{

}

//
void ScaffoldGroup::addLink(const ScaffoldLink& link)
{
    m_links.push_back(link);
}

// Calculate the longest overlap between any
// pair of nodes in the scaffold
// If two nodes are significantly overlapped, it
// may indicate that the root node is a repeat or connected
// to polymorphic nodes
int ScaffoldGroup::calculateLongestOverlap()
{
    int longestOverlap = 0;
    std::vector<ScaffoldLink>::iterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        // Compute interval for this link using closed coordionates [start, end]
        Interval interval_i(i->distance, i->getEndpoint() - 1);
        std::cout << "II: " << i->endpointID << " " << interval_i << "\n";

        std::vector<ScaffoldLink>::iterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            Interval interval_j(j->distance, j->getEndpoint() - 1);
            std::cout << "IJ: " << j->endpointID << " " << interval_j << "\n";
            if(Interval::isIntersecting(interval_i.start, interval_i.end,
                                        interval_j.start, interval_j.end)) {
                Interval intersection;
                Interval::intersect(interval_i.start, interval_i.end,
                                    interval_j.start, interval_j.end,
                                    intersection.start, intersection.end);

                int overlap = intersection.end - intersection.start + 1;
                if(overlap > longestOverlap)
                    longestOverlap = overlap;
            }
        }
    }

    return longestOverlap;
}

//
void ScaffoldGroup::computeBestOrdering()
{
    std::cout << "Compute ordering for" << m_pRootVertex->getID() << "\n";
    for(std::vector<ScaffoldLink>::iterator iter = m_links.begin();
                                            iter != m_links.end();
                                            ++iter)
    {
        std::cout << "\tlink: " << *iter << "\n";
    }
          
    // We compute the best ordering of the links for a given vertex
    // with a greedy algorithm as follows. Initially, all links
    // are unplaced. At each step, we test and score each unplaced
    // link as a candidate for being the next-closest contig
    // in the scaffold to the root. The lowest scoring contig
    // is removed from the unplaced list. This continues
    // until all contigs have been placed in the putative scaffold.
    LinkList unplacedLinks(m_links.begin(), m_links.end());
    std::vector<ScaffoldLink> placedLinks;

    int totalScore = 0;
    while(!unplacedLinks.empty())
    {
        int bestScore = std::numeric_limits<int>::max();
        LinkListIterator bestLink = unplacedLinks.end();

        for(LinkListIterator iter = unplacedLinks.begin();
                             iter != unplacedLinks.end();
                             ++iter)
        {
            int score = scoreLinkPlacement(*iter, unplacedLinks);
            std::cout << "Link score: " << score << "\n";

            if(score < bestScore)
            {
                bestScore = score;
                bestLink = iter;
            }
        }

        assert(bestLink != unplacedLinks.end());

        placedLinks.push_back(*bestLink);
        unplacedLinks.erase(bestLink);
        totalScore += bestScore;
    }

    std::cout << "Total ordering score: " << totalScore << "\n";
    std::cout << "Ordering: \n";
    for(std::vector<ScaffoldLink>::iterator iter = placedLinks.begin();
                                            iter != placedLinks.end();
                                            ++iter)
    {
        std::cout << "\tlink: " << *iter << "\n";
    }

    m_links = placedLinks;

    std::cout << "Longest overlap: " << calculateLongestOverlap() << "\n";
}

// Calculate the score of placing the given link before
// all the other links in the unplaced list
int ScaffoldGroup::scoreLinkPlacement(const ScaffoldLink& link,
                                      const LinkList& unplacedList)
{
    // Set the range of coordinates
    int end = link.getEndpoint();

    int score = 0;
    for(LinkList::const_iterator iter = unplacedList.begin();
                                 iter != unplacedList.end();
                                 ++iter)
    {
        if(link.endpointID != iter->endpointID)
        {
            // Compute the offset of this link based on input link preceding it
            // If the interval [start,end] includes the distance estimate
            // of the unplaced link, shift the estimate to the first valid
            // position after the placed contig (given by end - max overlap)
            int next_valid_position = end - m_maxOverlap - 1;
            if(iter->distance < next_valid_position)
            {
                score += next_valid_position - iter->distance;
            }
        }
    }
    return score;
}
