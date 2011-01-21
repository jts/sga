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
void ScaffoldGroup::addLink(const ScaffoldLink& link, ScaffoldVertex* pVertex)
{
    LinkVertexPair pair = {link, pVertex};
    m_links.push_back(pair);
}

// 
void ScaffoldGroup::resolveAmbiguity()
{
    double ambiguity_p = 0.01f;
    std::cout << "Checking " << m_pRootVertex->getID() << " for ambiguous links\n";
    LinkVectorIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        LinkVectorIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            bool isAmbiguous = areLinksAmbiguous(i->link, j->link, ambiguity_p);
            if(isAmbiguous)
            {
                double sum = i->pEndpoint->getEstCopyNumber() + j->pEndpoint->getEstCopyNumber();
                std::cout << "\tLinks " << i->link << " and " << j->link << " are ambiguous\n";
                std::cout << "\tECN I: "<< i->pEndpoint->getEstCopyNumber() << " ECN J: " << j->pEndpoint->getEstCopyNumber() << "\n";
                std::cout << "\tsum: " << sum << "\n";
            }
        }
    }
}

//
bool ScaffoldGroup::areLinksAmbiguous(const ScaffoldLink& linkA,
                                      const ScaffoldLink& linkB,
                                      double p)
{
    // Calculate the probabilty A comes before B and vice-versa
    double p_AB = calculateProbACloserThanB(linkA, linkB);
    double p_BA = 1.0f - p_AB;
    double best = std::max(p_AB, p_BA);
    double p_wrong = 1.0f - best;
    return p_wrong > p;
}



// Calculate the longest overlap between any
// pair of nodes in the scaffold
// If two nodes are significantly overlapped, it
// may indicate that the root node is a repeat or connected
// to polymorphic nodes
int ScaffoldGroup::calculateLongestOverlap()
{
    int longestOverlap = 0;
    LinkVectorIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        // Compute interval for this link using closed coordionates [start, end]
        Interval interval_i(i->link.distance, i->link.getEndpoint() - 1);
        std::cout << "II: " << i->link.endpointID << " " << interval_i << "\n";

        LinkVectorIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            Interval interval_j(j->link.distance, j->link.getEndpoint() - 1);
            std::cout << "IJ: " << j->link.endpointID << " " << interval_j << "\n";
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
    for(LinkVectorIterator iter = m_links.begin();
                                            iter != m_links.end();
                                            ++iter)
    {
        std::cout << "\tlink: " << iter->link << "\n";
    }
          
    // We compute the best ordering of the links for a given vertex
    // with a greedy algorithm as follows. Initially, all links
    // are unplaced. At each step, we test and score each unplaced
    // link as a candidate for being the next-closest contig
    // in the scaffold to the root. The lowest scoring contig
    // is removed from the unplaced list. This continues
    // until all contigs have been placed in the putative scaffold.
    LinkList unplacedLinks(m_links.begin(), m_links.end());
    LinkVector placedLinks;

    int totalScore = 0;
    while(!unplacedLinks.empty())
    {
        int bestScore = std::numeric_limits<int>::max();
        LinkListIterator bestLink = unplacedLinks.end();

        for(LinkListIterator iter = unplacedLinks.begin();
                             iter != unplacedLinks.end();
                             ++iter)
        {
            int score = scoreLinkPlacement(iter->link, unplacedLinks);
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
    for(LinkVectorIterator iter = placedLinks.begin(); iter != placedLinks.end(); ++iter)
    {
        std::cout << "\tlink: " << iter->link << "\n";
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
        if(link.endpointID != iter->link.endpointID)
        {
            std::cout << "P(" << link.endpointID << " < " << iter->link.endpointID << ") = " << calculateProbACloserThanB(link, iter->link) << "\n";
            // Compute the offset of this link based on input link preceding it
            // If the interval [start,end] includes the distance estimate
            // of the unplaced link, shift the estimate to the first valid
            // position after the placed contig (given by end - max overlap)
            int next_valid_position = end - m_maxOverlap - 1;
            if(iter->link.distance < next_valid_position)
            {
                score += next_valid_position - iter->link.distance;
            }
        }
    }
    return score;
}

// Given the links A and B, calculate the probabilty that A
// comes before B in the scaffold.
// In other words, if D_A and D_B are normal random variables, 
// what is the probability that D_B < D_A. The difference
// D_B - D_A is a random variable with mean u_A - u_B and variance
// v_A + v_B
double ScaffoldGroup::calculateProbACloserThanB(const ScaffoldLink& linkA,
                                                const ScaffoldLink& linkB)
{
    // calculate the mean and variance of the normal rv B - A
    double mean = linkA.distance - linkB.distance;
    double variance = pow(linkA.stdDev, 2.0f) + pow(linkB.stdDev, 2.0f);

    // The probabilty that A is closer than B is given by probabilty that
    // the random variable has a value less than zero.
    return normCDF(0.0f, mean, sqrt(variance));
}

double ScaffoldGroup::normCDF(double x, double mean, double sd)
{
    double t = (x - mean) / (sd * sqrt(2.0f));
    return 0.5f * (1 + erf(t)); 
}
