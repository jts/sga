#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include "SeqGraph.h"
#include "SeqVertex.h"

// Input objects
struct Contig
{
	VertexID id;
	int length;
	double coverage;
	std::string seq;

	friend std::istream& operator>> (std::istream& in, Contig& c)
	{
		in.ignore(1); // skip ">"
		in >> c.id >> c.length >> c.coverage; // read header
		in.ignore(1); // skip newline
		in >> c.seq; // read seq
		in.ignore(1); // place the ifstream at the next line
		return in;
	}
};


void loadContigVertices(SeqGraph& graph, int kmer, std::string filename)
{
	std::ifstream file(filename.c_str());
	assert(file.is_open());
	Contig c;
	while(file >> c)
	{
		SeqVertex* pSV = new SeqVertex(c.id, kmer, c.coverage, c.seq);
		graph.addVertex(pSV);
	}
}


int main(int argc, char** argv)
{
	(void)argv;
	(void)argc;
	
	SeqGraph sg;

	int argID = 1;
	int kmer = atoi(argv[argID++]);
	std::string contigFile(argv[argID++]);
	std::string adjFile(argv[argID++]);

	// Load contigs
	loadContigVertices(sg, kmer, contigFile);

	sg.stats();
}


void test()
{
	SeqGraph sg;
	Vertex* pV0 = new Vertex(0);
	Vertex* pV1 = new Vertex(1);
	Vertex* pV2 = new Vertex(2);
	Vertex* pV3 = new Vertex(3);
	Vertex* pV4 = new Vertex(4);
	sg.addVertex(pV0);
	sg.addVertex(pV1);
	sg.addVertex(pV2);
	sg.addVertex(pV3);
	sg.addVertex(pV4);

	sg.addEdge(pV0->getID(), pV2->getID(), ED_SENSE, EC_NATURAL);
	sg.addEdge(pV2->getID(), pV0->getID(), ED_ANTISENSE, EC_NATURAL);

	sg.addEdge(pV1->getID(), pV2->getID(), ED_SENSE, EC_NATURAL);
	sg.addEdge(pV2->getID(), pV1->getID(), ED_ANTISENSE, EC_NATURAL);

	sg.addEdge(pV2->getID(), pV3->getID(), ED_SENSE, EC_NATURAL);
	sg.addEdge(pV3->getID(), pV2->getID(), ED_ANTISENSE, EC_NATURAL);

	sg.addEdge(pV2->getID(), pV4->getID(), ED_SENSE, EC_NATURAL);
	sg.addEdge(pV4->getID(), pV2->getID(), ED_ANTISENSE, EC_NATURAL);

	sg.flip(4);
	sg.simplify();
	//sg.flip(4);

	sg.validate();
	sg.writeDot("blah");
}

