#include <stdio.h>
#include "SeqGraph.h"

int main(int argv, char** argc)
{
	(void)argv;
	(void)argc;
	
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

	sg.mergeVertices(0,2);
	sg.mergeVertices(1,2);

	sg.writeDot("blah");
}

