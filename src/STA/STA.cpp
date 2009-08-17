#include "STCommon.h"
#include "BWT.h"
#include "SuffixTree.h"
#include "ReadTable.h" 
#include "Timer.h"

void readDictTest()
{
	ReadTable rt;
	Read r1("r1", "AGATACAGAT");
	Read r2("r2", "GATACATAC");
	Read r3("r3", "ACATGATACAG");
	Read r4("r4", "TACAGATAT");

	rt.addRead(r1);
	rt.addRead(r2);
	rt.addRead(r3);
	rt.addRead(r4);

	SuffixArray sa1(0, r1.seq);
	SuffixArray sa2(1, r2.seq);
	SuffixArray sa3(2, r3.seq);
	SuffixArray sa4(3, r4.seq);
	
	SuffixArray sam1(sa1, sa2);
	SuffixArray sam2(sa3, sa4);
	SuffixArray sam3(sam1, sam2);

	sam1.validate(&rt);
	sam1.print(&rt);

	sam3.print(&rt);
	sam3.validate(&rt);
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
	SAVec* sav = new SAVec;
	Timer* t1 = new Timer("Building initial trees");
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		sav->push_back(SuffixArray(i, rt.getRead(i).seq));
		if(count % 10000 == 0)
			std::cout << "Built " << count << " trees\n";
		++count;
	}
	delete t1;

	Timer* t2 = new Timer("Merging trees");
	for(size_t j = 0; j < sav->size(); j+=2)
	{
		SuffixArray n((*sav)[j], (*sav)[j+1]);
	}
	delete t2;
	delete sav;
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


//
//
//
int main(int /*argc*/, char** argv)
{
	//simpleTest();
	(void)argv;
	//unitTest();
	//readDictTest();
	makeBWTFromReads(argv[1]);
	return 1;
}


