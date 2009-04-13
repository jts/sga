#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
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

struct AdjInfo
{
	EdgeVec edges;

	friend std::istream& operator>>(std::istream& in, AdjInfo& a)
	{
		VertexID root;
		std::string line;
		getline(in, line);

		// return if we've hit the end
		if(line == "")
			return in;
		
		std::stringstream parser(line);
		std::string bracket;
		std::string record;
		
		parser >> root; // parse the root
		parser >> bracket; // extract the opening bracket

		EdgeDir ed = ED_SENSE; // the first record is the sense record
		bool finished = false;

		while(!finished)
		{
			// Read the record
			parser >> record;

			if(record == "]")
			{
				if(ed == ED_SENSE)
				{
					ed = !ed;

					// extract the next bracket
					parser >> bracket;
				}
				else
				{
					finished = true;
				}
			}
			else
			{
				// Parse the record
				std::stringstream ss(record);
				
				VertexID id;
				int comp = -1;
				ss >> id;
				ss.ignore(1);
				ss >> comp;

				// Build the edge
				Edge e(root, id, ed, (EdgeComp)comp);
				a.edges.push_back(e);
			}
		}
		in.peek();
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

void loadContigEdges(SeqGraph& graph, std::string filename)
{
	std::ifstream file(filename.c_str());
	assert(file.is_open());
	AdjInfo a;
	while(file >> a)
	{
		EdgeVecIter iter = a.edges.begin();
		for(; iter != a.edges.end(); ++iter)
		{
			graph.addEdge(*iter);
		}
	}
	(void)graph;
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

	// Load verts and edges
	loadContigVertices(sg, kmer, contigFile);
	loadContigEdges(sg, adjFile);

	sg.stats();
	sg.validate();
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

