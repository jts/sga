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

// Construct a string from the scaffold
std::string ScaffoldRecord::generateString(const StringGraph* pGraph) const
{
    EdgeComp currComp = EC_SAME;

    // Starting from the root, join the sequence(s) of the scaffold
    // together along with the appropriate gaps/overlap
    Vertex* pVertex = pGraph->getVertex(m_rootID);
    assert(pVertex != NULL);
    
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

            // Calculate the strand this sequence is, on relative to the root
            if(link.getComp() == EC_REVERSE)
                currComp = !currComp;

            std::string toAppend = pVertex->getSeq().toString();
            if(currComp == EC_REVERSE)
            {
                toAppend = reverseComplement(toAppend);
            }

            if(reverseAll)
            {
                toAppend = reverse(toAppend);
            }


            std::string final;
            // Calculate the amount of overlap or gap characters to use
            if(link.distance < 0)
            {
                int remove = -1 * link.distance;
                final.append(10, 'N');
                final.append(toAppend.substr(remove));
            }
            else
            {
                final.append(link.distance, 'N');
                final.append(toAppend);
            }

            sequence.append(final);
        }

        if(reverseAll)
            sequence = reverse(sequence);
    }
    return sequence;
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
