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
std::string ScaffoldRecord::generateString(const StringGraph* pGraph, int minOverlap, int maxOverlap, double maxErrorRate) const
{
    int resolveMask = RESOLVE_GRAPH;// | RESOLVE_OVERLAP;


    // Starting from the root, join the sequence(s) of the scaffold
    // together along with the appropriate gaps/overlap
    Vertex* pVertex = pGraph->getVertex(m_rootID);
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
            const ScaffoldLink& link = m_links[i];

            pVertex = pGraph->getVertex(link.endpointID);

            // Calculate the strand this sequence is on relative to the root
            if(link.getComp() == EC_REVERSE)
                relativeComp = !relativeComp;

            //
            // Attempt to resolve the sequence of the link
            //
            std::string resolvedSequence;
            
            // Step 1, try to walk through the graph between the vertices
            bool resolved = false;
            if(resolveMask & RESOLVE_GRAPH)
            {
                resolved = graphResolve(pGraph, currID, link, resolvedSequence);
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
                    resolved = overlapResolve(sequence, toAppend, link, minOverlap, maxOverlap, maxErrorRate, resolvedSequence);

                // Step 3, just introduce a gap between the sequences
                if(!resolved)
                    introduceGap(toAppend, link, resolvedSequence);
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
                                  const ScaffoldLink& link, std::string& outExtensionString) const
{
    // Get the vertex to start the search from
    Vertex* pStartVertex = pGraph->getVertex(startID);
    Vertex* pEndVertex = pGraph->getVertex(link.endpointID);
    assert(pStartVertex != NULL && pEndVertex != NULL);

    SGWalkVector walks;
    SGSearch::findWalks(pStartVertex, pEndVertex, link.getDir(), 5000, 10000, walks);
    std::cout << "Found " << walks.size() << " paths that resolve the " << 
                  startID << " -> " <<link.endpointID << " link (de: " << link.distance << ")\n";

                 
    int NUM_STDDEV = 3;
    int threshold = NUM_STDDEV * link.stdDev;
    int numWalksConsistent = 0;
    std::string str;
    for(size_t i = 0; i < walks.size(); ++i)
    {
        int walkDistance = walks[i].getEndToStartDistance();
        int diff = abs(abs(link.distance - walkDistance));
        if(diff <= threshold)
        {
            ++numWalksConsistent;
            str = walks[i].getString(SGWT_EXTENSION);
            walks[i].print();
            std::cout << "Walk distance: " << walkDistance << " diff: " << diff << " threshold: " << threshold << "\n";
        }
    }
    std::cout << "Num walks satisfying distance: " << numWalksConsistent << "\n";

    if(numWalksConsistent == 1)
    {
        outExtensionString = str;
        return true;
    }
    else
    {
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
        upperBound = expectedOverlap + 3 * link.stdDev;
    else
        upperBound = maxOverlap;
    
    // Calculate the best match
    Match match;
    std::cout << "Searching for overlap of length: " << expectedOverlap << " max: " << upperBound << " sd: " << link.stdDev << "\n";
    bool overlapFound = OverlapTools::boundedOverlapDP(s1, s2, minOverlap, upperBound, maxErrorRate, match);
    if(overlapFound)
    {
        std::cout << "Overlap found, length: " << match.getMinOverlapLength() << "\n";
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
bool ScaffoldRecord::introduceGap(const std::string& contigString, const ScaffoldLink& link, std::string& out) const
{
    assert(out.empty());
    if(link.distance < 0)
    {
        // Truncate the string using the expected overlap and add a gap with a fixed number of Ns
        out.append(10, 'N');
        int expectedOverlap = -1 * link.distance;
        out.append(contigString.substr(expectedOverlap));
    }
    else
    {
        out.append(link.distance, 'N');
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
