//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSequenceCollection - A set of input
// sequences that are being scaffolded. Can
// either be held in a string graph or just a map
// from ID to sequence. Supports setting particular
// sequences as being visited
//
#include "ScaffoldSequenceCollection.h"

//
GraphSequenceCollection::GraphSequenceCollection(StringGraph* pGraph) : m_pGraph(pGraph)
{

}
        
// Returns the sequence with the given ID
std::string GraphSequenceCollection::getSequence(const std::string& id) const
{
    Vertex* pVertex = m_pGraph->getVertex(id);
    assert(pVertex != NULL);
    return pVertex->getSeq().toString();
}

// 
void GraphSequenceCollection::setPlaced(const std::string& id)
{
    Vertex* pVertex = m_pGraph->getVertex(id);
    assert(pVertex != NULL);
    pVertex->setColor(GC_BLACK);
}

// Small visitor function to write nodes of the graph that are not colored black
struct UnplacedVisitor
{
    UnplacedVisitor(std::ostream* pWriter, int minLength) : m_pWriter(pWriter), m_numUnplaced(0), m_minLength(minLength) {}

    void previsit(StringGraph*) {}
    
    bool visit(StringGraph*, Vertex* pVertex)
    {
        if(pVertex->getColor() == GC_BLACK)
            return false;
        
        if((int)pVertex->getSeqLen() >= m_minLength)
        {
            std::stringstream idss;
            idss << "unplaced-" << m_numUnplaced++;
            writeFastaRecord(m_pWriter, idss.str(), pVertex->getSeq().toString());
        }
        return false;
    }

    void postvisit(StringGraph*) {}

    std::ostream* m_pWriter;
    int m_numUnplaced;
    int m_minLength;
};

//
void GraphSequenceCollection::writeUnplaced(std::ostream* pWriter, int minLength)
{
    UnplacedVisitor uvisit(pWriter, minLength);
    m_pGraph->visit(uvisit);
}

