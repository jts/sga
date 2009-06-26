#include <iostream>
#include "Edge.h"
#include "Vertex.h"
#include "Bigraph.h"

int main(int argc, char** argv)
{
	(void)argc;
	(void)argv;

	typedef Edge<int> IntEdge;
	typedef Vertex<int, IntEdge> IntVertex;

	IntVertex* v1 = new IntVertex("1");
	IntVertex* v2 = new IntVertex("2");
	IntVertex* v3 = new IntVertex("3");

	typedef Bigraph<IntVertex> IntGraph;

	IntEdge e1("1", "2", ED_SENSE, EC_SAME, 10);
	IntEdge e2("2", "3", ED_SENSE, EC_SAME, 20);

	IntGraph g;
	g.addVertex(v1);
	g.addVertex(v2);
	g.addVertex(v3);
	v1->addEdge(e1);
	v2->addEdge(e2);
	g.writeDot("test.dot");
	return 0;
}


