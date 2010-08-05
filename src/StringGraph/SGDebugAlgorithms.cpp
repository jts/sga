//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGDebugAlgorithms - Methods used for the 
// development of algorithms.
//
#include "SGDebugAlgorithms.h"
#include "SGAlgorithms.h"
#include "SeqTrie.h"
#include "ErrorCorrect.h"
#include <algorithm>

//
// SGDebugEdgeClassificationVisitor - Collect statistics about the graph
// using debug information about simulated reads
//
void SGDebugEdgeClassificationVisitor::previsit(StringGraph* /*pGraph*/)
{
    num_good = 0;
    num_bad = 0;
    num_conflicted = 0;
    num_trusted = 0;
    num_nottrusted = 0;
}

//
bool SGDebugEdgeClassificationVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    assert(false && "Todo: reimplement");
    int curr_pos = 0; //pVertex->dbg_position;
    EdgePtrVec edges = pVertex->getEdges();

    if(pVertex->countEdges(ED_SENSE) > 1 || pVertex->countEdges(ED_ANTISENSE) > 1)
        num_conflicted++;

    for(size_t i = 0; i < edges.size(); ++i)
    {
        int edge_pos = 0; //edges[i]->getEnd()->dbg_position;
        int distance = abs(curr_pos - edge_pos);
        int overlap_len = edges[i]->getMatchLength();
        int sum = distance + overlap_len;
        
        if(sum == (int)pVertex->getSeq().length())
            ++num_good;
        else
            ++num_bad;

        if(edges[i]->isTrusted)
            num_trusted++;
        else
            num_nottrusted++;
    }
    return false;
}

// Remove all the marked edges
void SGDebugEdgeClassificationVisitor::postvisit(StringGraph* /*pGraph*/)
{
    printf("Num good: %d Num bad: %d Num conflicted: %d Num trusted: %d Num not trusted: %d\n", 
            num_good, num_bad, num_conflicted, num_trusted, num_nottrusted);
}

//
// SGDebugGraphCompareVisitor
// Compare the visited graph to the graph loaded
// in the constructor
//
SGDebugGraphCompareVisitor::SGDebugGraphCompareVisitor(std::string filename)
{
    m_showMissing = false;
    std::cout << "Loading graph " << filename << " for comparison\n";
    m_pCompareGraph = SGUtil::loadASQG(filename, 0, true);
}

//
SGDebugGraphCompareVisitor::~SGDebugGraphCompareVisitor()
{
    delete m_pCompareGraph;
    m_pCompareGraph = 0;
}

//
void SGDebugGraphCompareVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    m_numFound = 0;
    m_numMissing = 0;
    m_numMissingNull = 0;
    m_numWrong = 0;
    m_numContained = 0;
    m_numOpen = 0;
    m_numClosed = 0;
}

//
void SGDebugGraphCompareVisitor::postvisit(StringGraph*)
{
    std::cout << "[DebugSummary] NF: " << m_numFound << " NM: " << m_numMissing << 
                 " (" << m_numMissingNull << " NULL) " << " NW: " << m_numWrong <<
                 " NC: " << m_numContained << " NGO: " << m_numOpen << " NGC: " << m_numClosed << "\n";
    std::cout << "[DebugSummary] MissRate: " << double(m_numMissing) / double(m_numFound + m_numMissing) << "\n";
}


//
bool SGDebugGraphCompareVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    (void)pGraph;
    //compareTransitiveGroups(pGraph, pVertex);
    //compareInferredQuality(pGraph, pVertex);
    //compareOverlapQuality(pGraph, pVertex);
    //compareErrorRates(pGraph, pVertex);
    compareSplitGroups(pGraph, pVertex);
    /*
    if(!m_showMissing)
        summarize(pGraph, pVertex);
    else
        showMissing(pGraph, pVertex);
    */
    return false;
}

//
void SGDebugGraphCompareVisitor::showMissing(StringGraph* pGraph, Vertex* pVertex)
{
    (void)pGraph;
    // Retreive the vertex in the comparison graph
    Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
    if(pCompareVertex == NULL)
    {
        return;
    }

    EdgePtrVec compareEdges = pCompareVertex->getEdges();
    for(size_t i = 0; i < compareEdges.size(); ++i)
    {
        Edge* pCompareEdge = compareEdges[i];
        EdgeDesc ed = pCompareEdge->getDesc();
        if(!pVertex->hasEdge(ed))
        {
            std::cout << "MISSING!\t" << pCompareEdge->getMatchLength() << "\n";
            /*
            std::cout << "Edge: " << ed << " is missing in graph!";
            std::cout << "Vertex " << pVertex->getID() << " has edges:\n";
            EdgePtrVec actualEdges = pVertex->getEdges();
            for(size_t j = 0; j < actualEdges.size(); ++j)
            {
                std::cout << "\t" << *actualEdges[j] << "\n";
            }

            Vertex* pCompareEndpoint = pCompareEdge->getEnd();
            Vertex* pActualEndpoint = pGraph->getVertex(pCompareEndpoint->getID());
            if(pActualEndpoint == NULL)
                std::cout << "ENDPOINT IS NULL!\n";
            else
            {
                std::cout << "Endpoint " << pActualEndpoint->getID() << " has edges:\n";
                EdgePtrVec actualEdges = pActualEndpoint->getEdges();
                for(size_t j = 0; j < actualEdges.size(); ++j)
                {
                    std::cout << "\t" << *actualEdges[j] << "\n";
                }
            }
            */
        }
    }

}

//
void SGDebugGraphCompareVisitor::summarize(StringGraph* pGraph, Vertex* pVertex)
{
    // Retreive the vertex in the comparison graph
    Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
    if(pCompareVertex == NULL)
    {
        ++m_numContained;
        return;
    }

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];    
        TransitiveGroupCollection tgc = pVertex->computeTransitiveGroups(dir);
        if(tgc.numGroups() > 1)
            ++m_numOpen;
        else
            ++m_numClosed;
    }

    EdgePtrVec compareEdges = pCompareVertex->getEdges();
    for(size_t i = 0; i < compareEdges.size(); ++i)
    {
        Edge* pCompareEdge = compareEdges[i];
        EdgeDesc ed = pCompareEdge->getDesc();
        if(pVertex->hasEdge(ed))
        {
            ++m_numFound;
        }
        else
        {
            ++m_numMissing;
            Vertex* pEndpoint = pGraph->getVertex(pCompareEdge->getEndID());
            if(pEndpoint == NULL)
                ++m_numMissingNull;            
        }
    }

    EdgePtrVec actualEdges = pVertex->getEdges();
    for(size_t i = 0; i < actualEdges.size(); ++i)
    {
        Edge* pActualEdge = actualEdges[i];
        EdgeDesc ed = pActualEdge->getDesc();
        if(m_pCompareGraph->getVertex(pActualEdge->getEndID()) != NULL && !pCompareVertex->hasEdge(ed))
            ++m_numWrong;
    }

}

//
void SGDebugGraphCompareVisitor::compareSplitGroups(StringGraph* pGraph, Vertex* pVertex)
{
    // Retreive the vertex in the comparison graph
    Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
    if(pCompareVertex == NULL)
    {
        return;
    }

    MultiOverlap mo = pVertex->getMultiOverlap();    

    bool hasWrong = false;
    EdgePtrVec actualEdges = pVertex->getEdges();

    for(size_t i = 0; i < actualEdges.size(); ++i)
    {
        Edge* pActualEdge = actualEdges[i];
        EdgeDesc ed = pActualEdge->getDesc();
        if(m_pCompareGraph->getVertex(pActualEdge->getEndID()) != NULL && !pCompareVertex->hasEdge(ed))
        {
            hasWrong = true;
            break;
        }
    }

    std::string original = pVertex->getStr();
    std::string base = pCompareVertex->getStr();

    std::string d_str = getDiffString(base, original);
    std::cout << "ORG: " << original << "\n";
    std::cout << "CHT: " << base << "\n";
    std::cout << "DST: " << d_str << "\n";

    std::string consensus = ErrorCorrect::correctVertex(pGraph, pVertex, 3, 0.01);
    std::string c_diff_str = getDiffString(base, consensus);

    std::cout << "CNS: " << consensus << "\n";
    std::cout << "CDF: " << c_diff_str << "\n";

    //pVertex->setSeq(consensus);

    int numDiffsNC = countDifferences(original, base, base.size());
    int numDiffsEC = countDifferences(consensus, base, base.size());
    
    std::cout << "EC\t" << numDiffsNC << "\t" << numDiffsEC << "\t" << hasWrong << 
                 "\t" << pVertex->countEdges() << "\t" <<  "\n";

    if(numDiffsEC != 0)
    {
        for(size_t i = 0; i < consensus.size(); ++i)
        {
            if(consensus[i] != base[i])
                std::cout << "MM\t" << i << "\t" << numDiffsEC << "\t" << "\n";
        }
    }
}

void SGDebugGraphCompareVisitor::compareTransitiveGroups(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
    if(pCompareVertex == NULL)
    {
        return;
    }

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];    
        TransitiveGroupCollection actualTGC = pVertex->computeTransitiveGroups(dir);
        TransitiveGroupCollection compareTGC = pCompareVertex->computeTransitiveGroups(dir);

        if(actualTGC.numGroups() != compareTGC.numGroups())
        {
            printf("TGVMISMATCH\t%zu\t%zu\n", actualTGC.numGroups(), compareTGC.numGroups());
            for(size_t actualGroupIdx = 0; actualGroupIdx < actualTGC.numGroups(); ++actualGroupIdx)
            {
                Edge* pIrr = actualTGC.getGroup(actualGroupIdx).getIrreducible();
                size_t compareGroupIdx = compareTGC.findGroup(pIrr->getDesc());
                
                if(actualGroupIdx != compareGroupIdx)
                {
                    TransitiveGroup& compareGroup = compareTGC.getGroup(compareGroupIdx);
                    TransitiveGroup& actualGroup = actualTGC.getGroup(actualGroupIdx);

                    size_t actualGroupSize = actualGroup.numElements();
                    size_t compareGroupSize = compareGroup.numElements();

                    std::cout << "MPE\t" << *pIrr << "\t" << actualGroupIdx << "\t" << actualGroupSize << "\t"
                     << compareGroupIdx << "\t" << compareGroupSize << "\n";

                    /*
                    TransitiveGroup& compareGroup = compareTGC.getGroup(compareGroupIdx);
                    TransitiveGroup& actualGroup = actualTGC.getGroup(actualGroupIdx);

                    size_t actualGroupSize = actualGroup.numTransitive() + 1;
                    size_t compareGroupSize = compareGroup.numTransitive() + 1;

                    std::cout << "Compare collection:\n";
                    compareTGC.print();

                    std::cout << "Actual collection:\n";
                    actualTGC.print();

                    */
                    
                    // Compute the number of differences between the irreducible edge
                    // at the head of the incorrect group and all the elements of the 
                    // group it should be in
                    TransitiveGroup& expectedGroup = actualTGC.getGroup(compareGroupIdx);

                    MultiOverlap mo(pIrr->getEnd()->getID(), pIrr->getEnd()->getStr());

                    for(size_t i = 0; i < expectedGroup.numElements(); ++i)
                    {
                        Edge* pEdge = expectedGroup.getEdge(i);

                        // Set up the matches between the root vertex and i/j
                        // All coordinates are calculated from the point of view of pVertex
                        Match match_i = pIrr->getMatch();
                        Match match_j = pEdge->getMatch();

                        // Infer the match_ij based match_i and match_j
                        Match match_ij = Match::infer(match_i, match_j);
                        match_ij.expand();

                        // Convert the match to an overlap
                        Overlap ovr(pIrr->getEndID(), pEdge->getEndID(), match_ij);
                        mo.add(pEdge->getEnd()->getStr(), ovr);
                        int numDiff = ovr.match.countDifferences(pIrr->getEnd()->getStr(), pEdge->getEnd()->getStr());
                        std::cout << "MSE: " << pIrr->getEndID() << "\tInferred to: " << pEdge->getEndID() << " Num diff: ";
                        std::cout << numDiff << " error rate: " << double(numDiff) / double(ovr.match.getMinOverlapLength()) << "\n";
                    }
                }
            }
        }
    }
}

void SGDebugGraphCompareVisitor::compareErrorRates(StringGraph* pGraph, Vertex* pVertex)
{
    // Retreive the vertex in the comparison graph
    Vertex* pCompareVertex = m_pCompareGraph->getVertex(pVertex->getID());
    if(pCompareVertex == NULL)
    {
        std::cerr << "Compare vertex " << pVertex->getID() << " not found!\n";
        return;
    }
    EdgePtrVec compareEdges = pCompareVertex->getEdges();

    for(size_t i = 0; i < compareEdges.size(); ++i)
    {
        Edge* pCompareEdge = compareEdges[i];
        EdgeDesc ed = pCompareEdge->getDesc();
        bool is_missing = !pVertex->hasEdge(ed);

        //if(!pVertex->hasEdge(ed))
        {
            //std::cout << "Graph is missing edge " << ed << "\n";

            // Get the expected match between the sequences
            Match exp_match = pCompareEdge->getMatch();

            // Get the sequences of the endpoints for the current graph
            Vertex* pVertex2 = pGraph->getVertex(pCompareEdge->getEndID());

            std::string match1 = exp_match.coord[0].getSubstring(pVertex->getStr());
            std::string match2 = exp_match.coord[1].getSubstring(pVertex2->getStr());
            if(exp_match.isRC())
                match2 = reverseComplement(match2);
            
            std::string ms = is_missing ? "MISSING" : "FOUND";
            int num_diffs = countDifferences(match1, match2, match1.size());
            double mm_rate = double(num_diffs) / double(match1.size());
            (void)mm_rate;
            printf("%s\t%d\t%zu\t%lf\n", ms.c_str(), num_diffs, match1.size(), mm_rate);
            //std::cout << "SEQ1: " << pVertex->getStr() << "\n";
            //std::cout << "SEQ2: " << pVertex2->getStr() << "\n";
            //std::cout << "MATCH1: " << match1 << "\n";
            //std::cout << "MATCH2: " << match2 << "\n";
        }
    }
}
