#include <iostream>
#include "Edge.h"
#include "Vertex.h"
#include "Bigraph.h"

void graphTests();
void dnaStringTests();

int main(int argc, char** argv)
{
	(void)argc;
	(void)argv;
	dnaStringTests();
	return 0;
}

void dnaStringTests()
{
	std::string s1 = "GACACT";
	DNAString d1(s1);
	assert(s1 == d1.toString());
	std::cout << "S1: " << s1 << "\n";
	std::cout << "D1: " << d1.toString() << "\n";
	d1.reverse();
	std::cout << "Reversed: " << d1.toString() << "\n";

	d1 = "GGGACC";
	std::cout << "d1: " << d1.toString() << "\n";
	DNAString d2(d1);
	std::cout << "d2: " << d2.toString() << "\n"; 

	d1 = d2;
	std::cout << "d1: " << d1.toString() << "\n";
	d1 = d1;

	for(size_t i = 0; i < 2*d1.length(); ++i)
	{
		std::cout << i << " suffix: " << d1.getSuffix(i) << "\n";
		std::cout << i << " sufstr: " << d1.getSuffixString(i) << "\n";
	}
}

void graphTests()
{

	Vertex* v1 = new Vertex("1");
	Vertex* v2 = new Vertex("2");
	Vertex* v3 = new Vertex("3");

	Edge* pE1 = new Edge(v1, v2, ED_SENSE, EC_REVERSE);
	Edge* pE2 = new Edge(v2, v1, ED_SENSE, EC_REVERSE);
	Edge* pE3 = new Edge(v2, v3, ED_ANTISENSE, EC_SAME);
	Edge* pE4 = new Edge(v3, v2, ED_SENSE, EC_SAME);

	pE1->setTwin(pE2);
	pE2->setTwin(pE1);
	pE3->setTwin(pE4);
	pE4->setTwin(pE3);

	Bigraph g;
	g.addVertex(v1);
	g.addVertex(v2);
	g.addVertex(v3);
	g.addEdge(pE1);
	g.addEdge(pE2);
	g.addEdge(pE3);
	g.addEdge(pE4);
	g.validate();
	g.writeDot("before.dot");
	
	g.simplify();
	g.validate();
	g.writeDot("test.dot");
}


