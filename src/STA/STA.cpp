#include "STCommon.h"
#include "BWT.h"
#include "SuffixTree.h"
#include "ReadTable.h" 
#include "Timer.h"

void readDictTest()
{
	ReadTable rt;
	Read r1("r1", "ATA");
	Read r2("r2", "CTA");
	Read r3("r2", "ATA");

	rt.addRead(r1);
	rt.addRead(r3);
	rt.addRead(r2);
	SuffixArray sa;
	sa.initialize(rt);
	sa.sort(&rt);
	sa.validate(&rt);
}

void makeBWTFromReads(std::string file)
{
	std::ifstream in(file.c_str());
	ReadTable rt;

	std::string line;
	size_t count = 0;
	while(in >> line)
	{
		if(count % 2 == 1)
		{
			Read r("noname", line);
			rt.addRead(r);
		}

		if(count % 10000 == 0)
			std::cout << "Processed " << count << "\n";
		++count;
	}
	std::cout << "Loaded " << rt.getCount() << " reads\n";

	// Make initial suffix arrays
	SuffixArray sa;
	sa.initialize(rt);
	Timer* pT1 = new Timer("SuffixSort");
	sa.sort(&rt);
	delete pT1;
	sa.validate(&rt);
}


void makeBWTFromReads2(std::string file)
{
	std::ifstream in(file.c_str());
	ReadTable rt;

	std::string line;
	size_t count = 0;
	while(in >> line)
	{
		if(count % 2 == 1)
		{
			Read r("noname", line);
			rt.addRead(r);
			//std::cout << line << "\n";
			//reads.push_back(line);
			//BWT b(line);
			//b.printInfo();
		}

		if(count % 10000 == 0)
			std::cout << "Processed " << count << "\n";
		++count;
	}
	std::cout << "Loaded " << rt.getCount() << " reads\n";
	// Make initial suffix arrays
	typedef std::vector<SuffixArray> SAVec;

	count = 0;
	SAVec* inV = new SAVec;
	SAVec* outV = new SAVec;

	//
	//
	Timer* t1 = new Timer("Building initial trees");
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		inV->push_back(SuffixArray(i, rt.getRead(i).seq));
	}
	delete t1;

	//
	//
	Timer* t2 = new Timer("Merging trees");

	while(inV->size() > 1)
	{
		std::cout << "Merge: input has " << inV->size() << " elements\n";
		for(size_t j = 0; j < inV->size(); j+=2)
		{
			if(j != inV->size() - 1)
			{
				outV->push_back(SuffixArray((*inV)[j], (*inV)[j+1]));
			}
			else
			{
				outV->push_back((*inV)[j]);
			}
		}

		// Swap pointers
		inV->clear();
		SAVec* temp = inV;
		inV = outV;
		outV = temp;
	}
	std::cout << "Final vector has: " << inV->size() << " elements\n";
	inV->front().validate(&rt);
	inV->front().print(&rt);

	delete t2;
	delete inV;
	delete outV; 
	//BWT b(reads);
}

void unitTest()
{
	// Test construction of SAId's
	
	SAID said;
	uint32_t id1 = 3;
	uint32_t pos1 = 10;
	said.setID(id1);
	said.setPos(pos1);

	uint64_t id2 = 1249123132;
	uint8_t pos2 = 199;
	said.setPos(pos2);
	said.setID(id2);
	said.setID(id2);

	assert(said.getID() == id2);
	assert(said.getPos() == pos2);

}

