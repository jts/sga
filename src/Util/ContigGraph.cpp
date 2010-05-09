#include "ContigGraph.h"

//
// Build the vertices
//
void loadVertices(ContigGraph& graph, int /*kmer*/, std::string filename)
{
    std::ifstream file(filename.c_str());
    assert(file.is_open());
    Contig c;
    while(readCAF(file,c))
    {
        graph.addVertex(new ContigVertex(c.getID(), c));
    }
}

//
// Build the edges
//
void loadEdges(ContigGraph& graph, int overlap, std::string filename)
{
    std::ifstream file(filename.c_str());
    assert(file.is_open());
    AdjInfo a;
    while(file >> a)
    {
        graph.addEdge(new Edge(a.from, a.to, (EdgeDir)a.dir, (EdgeComp)a.comp));
    }
}

//
// Build the graph
//
ContigGraph* createContigGraph(int k, std::string contigsCAF, std::string adjCAF)
{
    ContigGraph* pGraph = new ContigGraph;
    loadVertices(*pGraph, k, contigsCAF);
    loadEdges(*pGraph, k - 1, adjCAF);
    return pGraph;
}

