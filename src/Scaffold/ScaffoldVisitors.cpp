//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVisitors - Functors that perform
// some operation on a ScaffoldVertex/Graph
//
#include "ScaffoldVisitors.h"
#include "ScaffoldRecord.h"
#include "ScaffoldGroup.h"
#include <limits.h>

//
// Returns true if the distance estimates for the two edges
// are close enough that the ordering cannot be resolved. Assumes
// that pXY is closer than pXZ
bool ScaffoldAlgorithms::areEdgesAmbiguous(ScaffoldEdge* pXY, ScaffoldEdge* pXZ)
{
    const double AMBIGUOUS_TOLERANCE = 2.0f;
    assert(pXY->getDistance() <= pXZ->getDistance());
    int translated = pXZ->getDistance() - static_cast<int>(AMBIGUOUS_TOLERANCE*pXZ->getStdDev());
    if(translated < pXY->getDistance())
        return true;
    else
        return false;
}

//
// Returns the edge properties inferred from the edges X->Y and X->Z
void ScaffoldAlgorithms::inferScaffoldEdgeYZ(ScaffoldEdge* pXY, ScaffoldEdge* pXZ,
                                             int& dist, double& sd, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp)
{
    comp = (pXY->getComp() == pXZ->getComp() ? EC_SAME : EC_REVERSE);
    dir_yz = !pXY->getTwin()->getDir(); // Opposite direction of X <- Y
    dir_zy = pXZ->getTwin()->getDir();  // Same direction as X <- Z since Z is after Y
    dist = pXZ->getDistance() - (pXY->getDistance() + pXY->getEnd()->getSeqLen());
    sd = sqrt(pow(pXY->getStdDev(), 2.0) + pow(pXZ->getStdDev(), 2.0));
}

// ScaffoldStatsVisitor
void ScaffoldStatsVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numVertices = 0;
    m_numEdges = 0;
}

bool ScaffoldStatsVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    ++m_numVertices;
    m_numEdges += pVertex->getNumEdges();

    return false;
}

void ScaffoldStatsVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    printf("Scaffold Stats -- Num vertices: %zu Num edges: %zu\n", m_numVertices, m_numEdges);
}

// ScaffoldAStatistic
ScaffoldAStatisticVisitor::ScaffoldAStatisticVisitor(double uniqueThreshold, 
                                                     double repeatThreshold) : m_uniqueThreshold(uniqueThreshold),
                                                                               m_repeatThreshold(repeatThreshold)
{
                                                                               
}

//
void ScaffoldAStatisticVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_sumUnique = 0;
    m_sumRepeat = 0;
    m_numUnique = 0;
    m_numRepeat = 0;
}

//
bool ScaffoldAStatisticVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() != SVC_REPEAT)
    {
        if(pVertex->getAStatistic() > m_uniqueThreshold)
        {
            pVertex->setClassification(SVC_UNIQUE);
            ++m_numUnique;
            m_sumUnique += pVertex->getSeqLen();
        }

        if(pVertex->getAStatistic() < m_repeatThreshold)
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numRepeat;
            m_sumRepeat += pVertex->getSeqLen();
        }
    }
    return false;   
}

//
void ScaffoldAStatisticVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    printf("A-statistic: classified %zu vertices as unique (%.2lf Mbp)\n", m_numUnique, (double)m_sumUnique / 1000000);
    printf("A-statistic: classified %zu vertices as repeat (%.2lf Mbp)\n", m_numRepeat, (double)m_sumRepeat / 1000000);
}

//
ScaffoldLinkValidator::ScaffoldLinkValidator(int maxOverlap, 
                                             double threshold) : m_maxOverlap(maxOverlap),
                                                                 m_threshold(threshold)
{

}

//
void ScaffoldLinkValidator::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numUnique = 0;
    m_numRepeat = 0;
}

//
bool ScaffoldLinkValidator::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            ScaffoldGroup group(pVertex, m_maxOverlap);
            for(size_t i = 0; i < edgeVec.size(); ++i)
            {
                group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());
            }

            group.resolveAmbiguity();
        }
    }
    return false;   
}

//
void ScaffoldLinkValidator::postvisit(ScaffoldGraph* /*pGraph*/)
{
    std::cerr << "Link validator done\n"; 
}

//
ScaffoldChainVisitor::ScaffoldChainVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{

}

//
void ScaffoldChainVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
}

// 
bool ScaffoldChainVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never try to make chains from a repeat
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            // Try to resolve this vertex edge set
            // by creating edges between the vertices that it
            // is linked to
            
            // Process the edges in order of length
            std::sort(edgeVec.begin(), edgeVec.end(), ScaffoldEdgePtrDistanceCompare);
            ScaffoldEdge* pXY = edgeVec[0];
            ScaffoldEdge* pXZ = edgeVec[1];
            assert(pXY->getDistance() <= pXZ->getDistance());
            ScaffoldVertex* pY = pXY->getEnd();
            ScaffoldVertex* pZ = pXZ->getEnd();

            if(pY->isRepeat() || pZ->isRepeat())
            {
                std::cerr << "Warning, skipping repeat\n";
                continue;
            }


            std::cout << "CV " << pVertex->getID() << "\t XY: " << *pXY << " XZ: " << *pXZ << "\n";

            // Check if the ordering of the contigs is ambiguous
            bool isAmbiguous = ScaffoldAlgorithms::areEdgesAmbiguous(pXY, pXZ);
            if(isAmbiguous)
            {
                std::cout << "\tEdges " << *pXY << " and " << *pXZ << 
                              " are ambiguously ordered\n";
            }

            // Infer the edge data for the edge Y->Z
            int dist;
            double sd;
            EdgeDir dir_yz;
            EdgeDir dir_zy; 
            EdgeComp comp;
            ScaffoldAlgorithms::inferScaffoldEdgeYZ(pXY, pXZ, dist, sd, dir_yz, dir_zy, comp);

            // Sanity checks
            bool isConsistent = (dist > -1*m_maxOverlap);
            if(!isConsistent)
            {
                std::cout << "\tEdge is not consistent with max overlap: " << dist << "\n";
            }

            // Check if an edge between Y->Z already exists
            ScaffoldEdge* pCheckEdge = pY->findEdgeTo(pZ->getID(), dir_yz, comp);
            if(pCheckEdge == NULL)
            {
                // Create the new edges
                ScaffoldLink linkYZ(pZ->getID(), dir_yz, comp, dist, sd, 0, pZ->getSeqLen(), SLT_INFERRED);
                ScaffoldLink linkZY(pY->getID(), dir_zy, comp, dist, sd, 0, pY->getSeqLen(), SLT_INFERRED);
                ScaffoldEdge* pYZ = new ScaffoldEdge(pZ, linkYZ);
                ScaffoldEdge* pZY = new ScaffoldEdge(pY, linkZY);
                pYZ->setTwin(pZY);
                pZY->setTwin(pYZ);

                std::cout << "\tCreating edge: " << *pYZ << " " << *pZY << "\n";
                pY->addEdge(pYZ);
                pZ->addEdge(pZY);
            }
            else
            {
                std::cout << "\tEdge already exists: " << *pCheckEdge << " inferred dist: " << dist << "\n";
            }

            // Remove the edges pXZ and pZX
            ScaffoldEdge* pZX = pXZ->getTwin();

            std::cout << "\tDeleting edge: " << *pXZ << "\n";
            pVertex->deleteEdge(pXZ); // deallocates
            pZ->deleteEdge(pZX); // deallocates
            
            pXZ = 0;
            pZX = 0;
            changed_graph = true;
        }
    }

    return changed_graph;
}

//
void ScaffoldChainVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{

}

//
bool ScaffoldMultiEdgeRemoveVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            pVertex->deleteEdgesAndTwins(dir);
        }
    }
    return false;
}

//
ScaffoldWriterVisitor::ScaffoldWriterVisitor(const std::string& filename)
{
    m_pWriter = createWriter(filename);
}

//
ScaffoldWriterVisitor::~ScaffoldWriterVisitor()
{
    delete m_pWriter;
}

//
void ScaffoldWriterVisitor::previsit(ScaffoldGraph* pGraph)
{
    pGraph->setVertexColors(GC_WHITE);
}

//
bool ScaffoldWriterVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    if(pVertex->getColor() == GC_RED)
        return false; //already output

    ScaffoldEdgePtrVector edges = pVertex->getEdges();
    
    if(edges.size() <= 1)
    {
        // Start of a chain found, traverse it to the other end
        size_t num_contigs = 0;
        size_t bases = 0; // number of bases in contigs
        size_t span = 0; // number of bases plus gaps

        pVertex->setColor(GC_RED);
    
        ScaffoldRecord record;
        record.setRoot(pVertex->getID());
        
        bases += pVertex->getSeqLen();
        num_contigs += 1;
        
        if(edges.size() == 1)
        {
            // Write the linked contigs
            ScaffoldEdge* pStartEdge = edges[0];
            ScaffoldEdge* pXY = pStartEdge;
            while(1)
            {
                record.addLink(pXY->getLink());
                ScaffoldVertex* pY = pXY->getEnd();
                if(pY->getColor() == GC_RED)
                    break; // loop found

                pY->setColor(GC_RED);
                bases += pY->getSeqLen();
                span += pY->getSeqLen() + pXY->getDistance();
                num_contigs += 1;

                // get the next direction
                EdgeDir nextDir = !pXY->getTwin()->getDir();
                ScaffoldEdgePtrVector nextEdges = pY->getEdges(nextDir);
                
                if(nextEdges.size() == 1)
                    pXY = nextEdges[0];
                else
                    break;
            }
        }
        record.writeScaf(m_pWriter);
        printf("Wrote scaffold with %zu components, %zu bases (%zu span)\n", num_contigs, bases, span);
    }
    return false;
}
