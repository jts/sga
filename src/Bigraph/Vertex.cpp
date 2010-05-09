//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#include "Vertex.h"
#include "Edge.h"
#include "TransitiveGroup.h"
#include <algorithm>

Vertex::~Vertex()
{
    EdgePtrVecIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        delete *iter;
        *iter = NULL;
    }

    if(m_pPairVertex != NULL)
        m_pPairVertex->clearPairVertex();
}

// Merging two string vertices has two parts
// First, the sequence of the vertex is extended
// by the the content of the edge label
// Then, all the edges that are pointing to this node
// must be updated to contain the extension of the vertex
void Vertex::merge(Edge* pEdge)
{
    Edge* pTwin = pEdge->getTwin();
    //std::cout << "Adding label to " << getID() << " str: " << pSE->getLabel() << "\n";

    // Merge the sequence
    std::string label = pEdge->getLabel();
    size_t label_len = label.size();
    pEdge->updateSeqLen(m_seq.length() + label.length());
    bool prepend = false;

    if(pEdge->getDir() == ED_SENSE)
    {
        m_seq.append(label);
    }
    else
    {
        label.insert(label.size(), m_seq);
        std::swap(m_seq, label);
        prepend = true;
    }

    pEdge->extendMatch(label_len);
    pTwin->extendMatchFullLength();

    // All the SeqCoords for the edges must have their seqlen field updated
    // Also, if we prepended sequence to this edge, all the matches in the 
    // SENSE direction must have their coordinates offset
    size_t newLen = m_seq.length();
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        Edge* pUpdateEdge = *iter;
        pUpdateEdge->updateSeqLen(newLen);
        if(prepend && pUpdateEdge->getDir() == ED_SENSE && pEdge != pUpdateEdge)
            pUpdateEdge->offsetMatch(label_len);
    }

    // Update the read count
    m_readCount += pEdge->getEnd()->getReadCount();

#ifdef VALIDATE
    VALIDATION_WARNING("Vertex::merge")
    validate();
#endif

}

void Vertex::validate() const
{
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        (*iter)->validate();
        /*
        std::string label = pSE->getLabel();
        StringVertex* pEnd = SV_CAST(pSE->getEnd());

        std::string vertSeq = pEnd->getSeq();
        EdgeDir suffixDir = (pSE->getComp() == EC_SAME) ? pSE->getDir() : !pSE->getDir();
        std::string vertSuffix = (suffixDir == ED_SENSE) ? vertSeq.substr(vertSeq.length() - label.length()) :
                                                                vertSeq.substr(0, label.length());

        if(pSE->getComp() == EC_REVERSE)
            vertSuffix = reverseComplement(vertSuffix);
        if(vertSuffix != label)
        {
            std::cerr << "Warning edge label " << label << " does not match vertex suffix " << vertSuffix << "\n";
        }
        */
    }
}

//
void Vertex::sortAdjListByID()
{
    EdgeIDComp comp;
    std::sort(m_edges.begin(), m_edges.end(), comp);
}

void Vertex::sortAdjListByLen()
{
    EdgeLenComp comp;
    std::sort(m_edges.begin(), m_edges.end(), comp);
}

// Ensure that this vertex has at most one edge per direction to any other vertex
void Vertex::makeUnique()
{
    // Sort the edge lists by length
    sortAdjListByLen();
    EdgePtrVec uniqueVec;
    makeUnique(ED_SENSE, uniqueVec);
    makeUnique(ED_ANTISENSE, uniqueVec);
    m_edges.swap(uniqueVec);
}

// Ensure each edge of the vertex is unique, assumes all vertices in the graph
// are initially colored white
void Vertex::makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec)
{
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        Edge* pEdge = *iter;
        if(pEdge->getDir() == dir)
        {
            Vertex* pY = pEdge->getEnd();
            if(pY->getColor() == GC_BLACK)
            {
                std::cerr << getID() << " has a duplicate edge to " << pEdge->getEndID() << " in direction " << dir << "\n";
                // This vertex is the endpoint of some other (longer) edge
                // Delete the edge and remove it from the twin
                Edge* pTwin = pEdge->getTwin();
                Vertex* pPartner = pEdge->getEnd();
                pPartner->removeEdge(pTwin);
                delete pEdge;
                pEdge = NULL;
                delete pTwin;
                pTwin = NULL;
            }
            else
            {
                assert(pY->getColor() == GC_WHITE);
                pY->setColor(GC_BLACK);
                uniqueVec.push_back(*iter);
            }
        }
    }

    for(EdgePtrVecIter iter = uniqueVec.begin(); iter != uniqueVec.end(); ++iter)
        (*iter)->getEnd()->setColor(GC_WHITE);
}


// Compute the transitive groups for this vertex
// A transitive group is a set of edges s.t. one
// edge is irreducible and the other edges in the set
// are transitive w.r.t. the irreducible edge.
TransitiveGroupCollection Vertex::computeTransitiveGroups(EdgeDir dir)
{
    WARN_ONCE("Compute transitive groups could be better");
    TransitiveGroupCollection tgc(this, dir);

    // The edges must be sorted by the length of the overlap
    sortAdjListByLen();

    static const size_t FUZZ = 10; // see myers
    
    EdgePtrVec edges = getEdges(dir);

    if(edges.size() == 0)
        return tgc;

    // Set all the vertices connected to this as unvisited
    for(size_t i = 0; i < edges.size(); ++i)
        (edges[i])->getEnd()->setColor(GC_GRAY);

    Edge* pLongestEdge = edges.back();
    size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;
    
    // Stage 1
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pVWEdge = edges[i];
        Vertex* pWVert = pVWEdge->getEnd();

        EdgeDir transDir = !pVWEdge->getTwinDir();
        if(pWVert->getColor() == GC_GRAY)
        {
            EdgePtrVec w_edges = pWVert->getEdges(transDir);
            for(size_t j = 0; j < w_edges.size(); ++j)
            {
                Edge* pWXEdge = w_edges[j];
                size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
                if(trans_len <= longestLen)
                {
                    if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                    {
                        // X is the endpoint of an edge of V, therefore it is transitive
                        pWXEdge->getEnd()->setColor(GC_BLACK);
                    }
                }
                else
                {
                    break;
                }
            }
        }
    }
    
    /*
    // Stage 2
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pVWEdge = edges[i];
        Vertex* pWVert = pVWEdge->getEnd();

        //std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
        //std::cout << pWVert->getID() << " w_edges: \n";
        EdgeDir transDir = !pVWEdge->getTwinDir();
        EdgePtrVec w_edges = pWVert->getEdges(transDir);
        for(size_t j = 0; j < w_edges.size(); ++j)
        {
            //std::cout << "    edge: " << *w_edges[j] << "\n";
            Edge* pWXEdge = w_edges[j];
            size_t len = pWXEdge->getSeqLen();

            if(len < FUZZ || j == 0)
            {
                if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                {
                    // X is the endpoint of an edge of V, therefore it is transitive
                    pWXEdge->getEnd()->setColor(GC_BLACK);
                    ++marked_verts;
                    //std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << " in stage 2\n";
                }
            }
            else
                break;
        }
    }
    */
    
    // At this point the vertices of transitive edges are marked black
    // and the vertices of irredicuble edges are GRAY. Create the transitive groups
    for(size_t i = 0; i < edges.size(); ++i)
    {
        if(edges[i]->getEnd()->getColor() == GC_GRAY)
        {
            Edge* pIrreducibleEdge = edges[i];
            tgc.createGroup(pIrreducibleEdge);
        }
        
        // Reset colors
        edges[i]->getEnd()->setColor(GC_WHITE);
    }

    // Fill out the newly created group
    for(size_t i = 0; i < tgc.numGroups(); ++i)
    {
        TransitiveGroup& group = tgc.getGroup(i);
        Edge* pIrreducibleEdge = group.getIrreducible();
        Vertex* pIrreducibleVert = pIrreducibleEdge->getEnd();
        EdgeDir transDir = !pIrreducibleEdge->getTwinDir();
        EdgePtrVec irrEdges = pIrreducibleVert->getEdges(transDir);
        
        // Mark the reachable edges as black
        for(size_t j = 0; j < irrEdges.size(); ++j)
        {
            irrEdges[j]->getEnd()->setColor(GC_BLACK);
        }
    
        // Add the black edges reachable from *this to the group
        for(size_t j = 0; j < edges.size(); ++j)
        {
            if(edges[j]->getEnd()->getColor() == GC_BLACK)
            {
                group.add(edges[j]);
            }
        }

        for(size_t j = 0; j < irrEdges.size(); ++j)
        {
            irrEdges[j]->getEnd()->setColor(GC_WHITE);
        }
    }    
    return tgc;
}

// Get a multioverlap object representing the overlaps for this vertex
MultiOverlap Vertex::getMultiOverlap() const
{
    MultiOverlap mo(getID(), getSeq());
    for(size_t i = 0; i < m_edges.size(); ++i)
    {
        Edge* pEdge = m_edges[i];
        mo.add(pEdge->getEnd()->getSeq(), pEdge->getOverlap());
    }
    return mo;
}

// Get a SeqTrie of the overlaps of this vertex
void Vertex::fillTries(double p_error, SeqTrie* pSenseTrie, SeqTrie* pAntisenseTrie) const
{
    double lp = log(p_error);

    if(m_edges.empty())
        return;
    
    for(size_t i = 0; i < m_edges.size(); ++i)
    {
        Edge* pEdge = m_edges[i];
        std::string overlapped = pEdge->getTwin()->getMatchStr();
        if(pEdge->getComp() == EC_REVERSE)
            overlapped = reverseComplement(overlapped);

        if(pEdge->getDir() == ED_SENSE)
        {
            overlapped = reverse(overlapped);
            pSenseTrie->insert(overlapped, lp);
        }
        else
        {
            pAntisenseTrie->insert(overlapped, lp);
        }
    }
}


// Get the inferred probabilty of each base for this sequence
// base on the multioverlap
QualityVector Vertex::getInferredQuality() const
{
    MultiOverlap mo = getMultiOverlap();
    QualityVector out;
    for(size_t i = 0; i < m_seq.size(); ++i)
    {
        DNADouble ap = mo.calcAlphaProb(i);
        out.add(ap);
    }
    return out;
}

//
std::string Vertex::getInferredConsensus() const
{
    WARN_ONCE("Vertex consensus is heurestic"); 
    MultiOverlap mo = getMultiOverlap();
    std::string out;
    out.reserve(m_seq.size());
    for(size_t i = 0; i < m_seq.size(); ++i)
    {
        AlphaCount ac = mo.calcAlphaCount(i);
        char b = ac.getMaxBase();
        if(ac.get(b) > ac.get(m_seq[i]))
            out.push_back(b);
        else
            out.push_back(m_seq[i]);
    }
    return out;

    /*
    MultiOverlap mo = getMultiOverlap();
    std::string out;
    out.reserve(m_seq.size());
    for(size_t i = 0; i < m_seq.size(); ++i)
    {
        AlphaProb ap = mo.calcAlphaProb(i);
        char b = ap.getMaxBase();
        out.push_back(b);
    }
    return out;
    */
}


// Get the prior probability of each base of the sequence
// based on either quality scores or a flat prior
QualityVector Vertex::getPriorQuality() const
{
    WARN_ONCE("Using flat prior of 0.99");
    QualityVector out;
    for(size_t i = 0; i < m_seq.size(); ++i)
    {
        DNADouble ap;
        double p_correct_match = 0.99;
        double p_correct_mm = (1.0f - p_correct_match) / 3;
        for(size_t j = 0; j < DNA_ALPHABET_SIZE; ++j)
        {
            char b = ALPHABET[j];
            if(b == m_seq[i])
            {
                ap.set(b, log(p_correct_match));
            }
            else
            {
                ap.set(b, log(p_correct_mm));
            }
        }
        out.add(ap);
    }
    return out;
}

// Add an edge
void Vertex::addEdge(Edge* ep)
{
    assert(ep->getStart() == this);

#ifdef VALIDATE
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == ep->getEndID())
        {
            std::cout << "Attempted to add duplicate edge with ID: " << ep->getEndID() 
                        << " to vertex: " << ep->getStartID() << "\n";
            std::cout << "Added in desc: " << ep->getDesc() << " curr desc: " << (*iter)->getDesc() << "\n";
            assert(false);
        }
    }
#endif
    m_edges.push_back(ep);
}

// Remove an edge from the edge list of the vertex
// and free the storage. This does not remove the twin
void Vertex::deleteEdge(Edge* pEdge)
{
    removeEdge(pEdge);
    delete pEdge;
    pEdge = NULL;
}

// Remove an edge but do not destroy it
void Vertex::removeEdge(Edge* pEdge)
{
    EdgePtrVecIter iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if(*iter == pEdge)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);    
}

//
void Vertex::removeEdge(const EdgeDesc& ed)
{
    EdgePtrVecIter iter = findEdge(ed);
    if(iter == m_edges.end())
    {
        std::cout << "EDGE NOT FOUND: " << ed << "\n";
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);
}

// Delete all the edges, and their twins, from this vertex
void Vertex::deleteEdges()
{
    EdgePtrVecIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        Edge* pEdge = *iter;
        Edge* pTwin = pEdge->getTwin();
        Vertex* pPartner = pEdge->getEnd();
        pPartner->removeEdge(pTwin);
        delete pEdge;
        pEdge = NULL;
        delete pTwin;
        pTwin = NULL;
        *iter = NULL;
    }
    m_edges.clear();
}

// Delete edges that are marked
// This only deletes the edge and not its twin
void Vertex::sweepEdges(GraphColor c)
{
    EdgePtrVecIter iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        Edge* pEdge = *iter;
        if(pEdge->getColor() == c)
        {
            delete pEdge;
            pEdge = NULL;
            iter = m_edges.erase(iter);
        }
        else
            ++iter;
    }
}

// Return the iterator to the edge matching edgedesc
EdgePtrVecIter Vertex::findEdge(const EdgeDesc& ed)
{
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDesc() == ed)
            return iter;
    }
    return m_edges.end();
}

//
EdgePtrVecConstIter Vertex::findEdge(const EdgeDesc& ed) const
{
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDesc() == ed)
            return iter;
    }
    return m_edges.end();
}

// Check for the presence of an edge
bool Vertex::hasEdge(Edge* pEdge) const
{
    return hasEdge(pEdge->getDesc());
}

//
bool Vertex::hasEdge(const EdgeDesc& ed) const
{
    return findEdge(ed) != m_edges.end();
}

//
bool Vertex::hasEdgeTo(const Vertex* pY) const
{
    assert(pY != NULL);
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        if((*iter)->getEnd() == pY)
            return true;
    return false;
}

// Return the edge matching the descriptions
Edge* Vertex::getEdge(const EdgeDesc& ed)
{
     EdgePtrVecIter i = findEdge(ed);
     assert(i != m_edges.end());
     return *i;
}

// Find edges to the specified vertex
EdgePtrVec Vertex::findEdgesTo(VertexID id)
{
    EdgePtrVecConstIter iter = m_edges.begin();
    EdgePtrVec outEdges;
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id)
            outEdges.push_back(*iter);
    }
    return outEdges;
}

//
void Vertex::setPairVertex(Vertex* pPair) 
{ 
    m_pPairVertex = pPair; 
}

//
void Vertex::clearPairVertex() 
{ 
    m_pPairVertex = NULL; 
}


//
// Get the edges in a particular direction
// This preserves the ordering of the edges
//
EdgePtrVec Vertex::getEdges(EdgeDir dir) const
{
    EdgePtrVecConstIter iter = m_edges.begin();
    EdgePtrVec outEdges;
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() == dir)
            outEdges.push_back(*iter);
    }
    return outEdges;
}


// Get the edges
EdgePtrVec Vertex::getEdges() const
{
    EdgePtrVec outEdges(m_edges.begin(), m_edges.end());
    return outEdges;    
}

void Vertex::setEdgeColors(GraphColor c) 
{ 
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        (*iter)->setColor(c);
}

// Count the edges
// This function is not necessarily constant time
size_t Vertex::countEdges() const
{ 
    return m_edges.size(); 
}

//
size_t Vertex::countEdges(EdgeDir dir)
{
    EdgePtrVec ev = getEdges(dir);
    return ev.size();
}

// Return the amount of memory this vertex is using, in bytes
size_t Vertex::getMemSize() const
{
    return sizeof(*this) + (m_edges.size() * sizeof(Edge*) + (m_seq.length() * sizeof(char)));
}


// Output edges in graphviz format
void Vertex::writeEdges(std::ostream& out, int dotFlags) const
{
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        if(dotFlags & DF_UNDIRECTED)
        {
            if((*iter)->getStartID() < (*iter)->getEndID())
            {
                out << "\"" << (*iter)->getStart() << "\" -- \"" << (*iter)->getEnd() << "\"";
            }
        }
        else
        {
            out << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
            std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
            std::string label = ((*iter)->getComp() == EC_SAME) ? "S" : "F";
            out << "\" [color=\"" << color << "\" ";
            out << "label=\"" << (*iter)->getLabel() << "\"];";
        }
        out << "\n";
    }
}

bool EdgeIDComp::operator()(const Edge* pA, const Edge* pB) 
{
       return pA->getEndID() < pB->getEndID();
}

// Compare string edge points by length
bool EdgeLenComp::operator()(const Edge* pA, const Edge* pB)
{
    return pA->getSeqLen() < pB->getSeqLen();
}

