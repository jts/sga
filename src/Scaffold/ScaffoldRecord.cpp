//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldRecord - A scaffold consisting of a 
// starting component and a vector of ordered
// links
//
#include "ScaffoldRecord.h"
#include "OverlapTools.h"
#include "SGSearch.h"

//#define DEBUGRESOLVE 1

//
ScaffoldRecord::ScaffoldRecord() 
{

}

//
void ScaffoldRecord::setRoot(const std::string& root)
{
    m_rootID = root;
}

//
void ScaffoldRecord::addLink(const ScaffoldLink& link)
{
    m_links.push_back(link);
}

//
size_t ScaffoldRecord::getNumComponents() const
{
    if(m_rootID.empty())
        return 0;
    
    return 1 + m_links.size();
}

// Construct a string from the scaffold
std::string ScaffoldRecord::generateString(const StringGraph* pGraph, int minOverlap, 
                                           int maxOverlap, double maxErrorRate, int resolveMask, 
                                           int minGapLength, ResolveStats* pStats) const
{
    pStats->numScaffolds += 1;

    // Starting from the root, join the sequence(s) of the scaffold
    // together along with the appropriate gaps/overlap
    Vertex* pVertex = pGraph->getVertex(m_rootID);
    pVertex->setColor(GC_BLACK);

    assert(pVertex != NULL);
    
    EdgeComp relativeComp = EC_SAME;
    EdgeComp prevComp = EC_SAME;

    std::string currID = m_rootID;
    std::string sequence = pVertex->getSeq().toString();

    // Add in the sequences of linked contigs, if any
    if(m_links.size() > 0)
    {
        EdgeDir rootDir = m_links[0].getDir();
        
        // If this scaffold grows in the antisense direction,
        // we reverse every component and perform appends of the reversed
        // parts. After the scaffold is constructed we reverse again
        // to obtain the final scaffold in the desired orientation
        bool reverseAll = (rootDir == ED_ANTISENSE);
        if(reverseAll)
            sequence = reverse(sequence);

        for(size_t i = 0; i < m_links.size(); ++i)
        {   
            pStats->numGapsAttempted += 1;
            const ScaffoldLink& link = m_links[i];

            pVertex = pGraph->getVertex(link.endpointID);
            pVertex->setColor(GC_BLACK);

            // Calculate the strand this sequence is on relative to the root
            if(link.getComp() == EC_REVERSE)
                relativeComp = !relativeComp;

            //
            // Attempt to resolve the sequence of the link
            //
            std::string resolvedSequence;
            
            // Step 1, try to walk through the graph between the vertices
            bool resolved = false;
            if(resolveMask & RESOLVE_GRAPH_BEST || resolveMask & RESOLVE_GRAPH_UNIQUE)
            {
                resolved = graphResolve(pGraph, currID, link, resolveMask, pStats, resolvedSequence);
                /*
                std::cout << "GR " << link.endpointID << " PC: " << prevComp << " LC: " << 
                             link.getComp() << " RC: " << relativeComp << " LD: " << 
                             link.getDir() << " resolved: " << resolved << "\n";
                */
                
                if(resolved)
                {
                    // The returned sequence is wrt currID. If we flipped the sequence of currID, we must
                    // flip this sequence
                    if(prevComp == EC_REVERSE)
                        resolvedSequence = reverseComplement(resolvedSequence);

                    if(reverseAll)
                        resolvedSequence = reverse(resolvedSequence);
                    pStats->numGapsResolved += 1;
                }
            }

            // Step 2, try to resolve a predicted overlap between current sequence
            // and the sequence of the linked contig
            if(!resolved)
            {
                // Get the sequence that should be potentially appended in
                std::string toAppend = pVertex->getSeq().toString();
                if(relativeComp == EC_REVERSE)
                    toAppend = reverseComplement(toAppend);
                if(reverseAll)
                    toAppend = reverse(toAppend);
                
                //
                if(!resolved && link.distance < 0 && resolveMask & RESOLVE_OVERLAP)
                {
                    resolved = overlapResolve(sequence, toAppend, link, minOverlap, maxOverlap, maxErrorRate, resolvedSequence);
                    if(resolved)
                    {
                        pStats->numGapsResolved += 1;
                        pStats->overlapFound += 1;
                    }
                    else
                    {
                        pStats->overlapFailed += 1;
                    }
                }

                // Step 3, just introduce a gap between the sequences
                if(!resolved)
                    introduceGap(minGapLength, toAppend, link, resolvedSequence);

            }

            sequence.append(resolvedSequence);
            currID = link.endpointID;
            prevComp = relativeComp;
        }

        if(reverseAll)
            sequence = reverse(sequence);
    }
    return sequence;
}

// Attempt to resolve a scaffold link by finding a walk through the graph linking the two vertices
bool ScaffoldRecord::graphResolve(const StringGraph* pGraph, const std::string& startID, 
                                  const ScaffoldLink& link, int resolveMask, 
                                  ResolveStats* pStats, std::string& outExtensionString) const
{
    // Get the vertex to start the search from
    Vertex* pStartVertex = pGraph->getVertex(startID);
    Vertex* pEndVertex = pGraph->getVertex(link.endpointID);
    assert(pStartVertex != NULL && pEndVertex != NULL);

    double NUM_STDDEV = 3.0f;
    int threshold = static_cast<int>(NUM_STDDEV * link.stdDev);
    int maxDistance = link.distance + threshold;
    int maxExtensionDistance = maxDistance + pEndVertex->getSeqLen();
    SGWalkVector walks;
    SGSearch::findWalks(pStartVertex, pEndVertex, link.getDir(), maxExtensionDistance, 10000, true, walks);

    int numWalksValid = 0;
    int numWalksClosest = 0;
    int selectedIdx = -1;
    int closestDist = std::numeric_limits<int>::max();
    
    // Select the closest walk to the distance estimate
    for(size_t i = 0; i < walks.size(); ++i)
    {
        int walkDistance = walks[i].getEndToStartDistance();
        int diff = abs(abs(link.distance - walkDistance));
        if(diff <= threshold)
        {

#ifdef DEBUGRESOLVE
            std::cout << "Walk distance: " << walkDistance << " diff: " << diff << " threshold: " << threshold << " close: " << closestDist << "\n";
#endif
            ++numWalksValid;
            if(diff < closestDist)
            {
                selectedIdx = i;
                closestDist = diff;
                numWalksClosest = 1;
            }
            else if(diff == closestDist)
            {
                numWalksClosest += 1;
            }

        }
    }

    // Choose the best path, if any, depending on the algorithm to use
    bool useWalk = false;

    if(numWalksValid > 0)
    {
        if(resolveMask & RESOLVE_GRAPH_BEST)
        {
            // If the unique flag is not set, or we only have 1 closest walk, select it
            if(!(resolveMask & RESOLVE_GRAPH_UNIQUE) || numWalksClosest == 1)
                useWalk = true;
            else if((resolveMask & RESOLVE_GRAPH_UNIQUE) && numWalksClosest > 1)
                pStats->graphWalkTooMany += 1;
        }
        else
        {
            if(numWalksValid == 1)
                useWalk = true;
            else if(numWalksValid > 1)
                pStats->graphWalkTooMany += 1;
        }
    }

#ifdef DEBUGRESOLVE    
    std::cout << "Num walks: " << walks.size() << " Num valid: " << numWalksValid << " Num closest: " << numWalksClosest << " using: " << useWalk << "\n";
#endif

    if(useWalk)
    {
        assert(selectedIdx != -1);
        outExtensionString = walks[selectedIdx].getString(SGWT_EXTENSION);
        pStats->graphWalkFound += 1;

        // Mark all vertices in the walk as visited
        VertexPtrVec vertexPtrVector = walks[selectedIdx].getVertices();
        for(size_t i = 0; i < vertexPtrVector.size(); ++i)
            vertexPtrVector[i]->setColor(GC_BLACK);
        return true;
    }
    else
    {
        if(numWalksValid == 0)
            pStats->graphWalkNoPath += 1;
        assert(outExtensionString.empty());
        return false;
    }
}

// Attempt to resolve a predicted overlap between s1 and s2
// Returns true if there overlap was found and the overhang of s2 is placed in outString
bool ScaffoldRecord::overlapResolve(const std::string& s1, const std::string& s2, 
                                    const ScaffoldLink& link, int minOverlap, int maxOverlap, 
                                    double maxErrorRate, std::string& outString) const
{
    // Attempt to find an overlap between these sequences
    int expectedOverlap = -1 * link.distance;

    // If the maximum overlap was not set, set it to the expected overlap * 3 stddev
    int upperBound = 0;
    if(maxOverlap == -1)
        upperBound = static_cast<int>(expectedOverlap + 3.0f * link.stdDev);
    else
        upperBound = maxOverlap;
    
    // Calculate the best match
    Match match;
    bool overlapFound = OverlapTools::boundedOverlapDP(s1, s2, minOverlap, upperBound, maxErrorRate, match);
    if(overlapFound)
    {
        SeqCoord overlapCoord = match.coord[1];
        SeqCoord overhangCoord = overlapCoord.complement();
        outString = overhangCoord.getSubstring(s2);
        return true;
    }
    else
    {
        return false;
    }
}

// Resolve a link with a gap
bool ScaffoldRecord::introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& out) const
{
    assert(out.empty());
    if(link.distance < 0)
    {
        // Truncate the string using the expected overlap and add a gap with a fixed number of Ns
        out.append(minGapLength, 'N');
        int expectedOverlap = -1 * link.distance;
        assert(expectedOverlap < (int)contigString.length());
        out.append(contigString.substr(expectedOverlap));
    }
    else
    {
        int gap = std::max(link.distance, minGapLength);
        out.append(gap, 'N');
        out.append(contigString);
    }
    return true;
}

//
void ScaffoldRecord::parse(const std::string& text)
{
    StringVector fields = split(text, '\t');
    assert(fields.size() >= 1);

    m_rootID = fields[0];
    for(size_t i = 1; i < fields.size(); ++i)
    {
        ScaffoldLink link;
        link.parse(fields[i]);
        m_links.push_back(link);
    }
}

//
void ScaffoldRecord::writeScaf(std::ostream* pWriter)
{
    *pWriter << m_rootID;
    for(size_t i = 0; i < m_links.size(); ++i)
        *pWriter << "\t" << m_links[i];
    *pWriter << "\n";
}
