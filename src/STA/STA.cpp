#include "STCommon.h"
#include "BWT.h"
#include "SuffixTree.h"
#include "ReadTable.h" 

void readDictTest()
{
	ReadTable rt;
	Read r1("r1", "AGATACAGAT");
	Read r2("r2", "GATACATAC");

	rt.addRead(r1);
	rt.addRead(r2);

	SuffixArray sa1(0, r1.seq);
	SuffixArray sa2(1, r2.seq);
	SuffixArray sam(sa1, sa2);
	sam.validate(&rt);
	sam.print(&rt);
}

void makeBWTFromReads(std::string file)
{
	StringVector reads;
	std::ifstream in(file.c_str());
	std::string line;
	size_t count = 0;
	while(in >> line)
	{
		if(count % 2 == 1)
		{
			//std::cout << line << "\n";
			//reads.push_back(line);
			//BWT b(line);
			//b.printInfo();
			return;
		}

		if(count % 10000 == 0)
			std::cout << "Processed " << count << "\n";
		++count;
	}
	std::cout << "Loaded " << reads.size() << " reads\n";
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
	readDictTest();
	//makeBWTFromReads(argv[1]);
	return 1;
}


