#include <stdio.h>
#include "SeqGraph.h"

int main(int argv, char** argc)
{
	(void)argv;
	(void)argc;
	
	SeqGraph sg;

	IVertex* pV0 = new IVertex(0);
	IVertex* pV1 = new IVertex(1);
	IVertex* pV2 = new IVertex(2);
	IVertex* pV3 = new IVertex(3);
	IVertex* pV4 = new IVertex(4);
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


	sg.writeDot("blah");
}

