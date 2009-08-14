#include "STCommon.h"
#include "BWT.h"
#include "SuffixTree.h"

//
//
//
void dnaTest()
{
	std::string ref = "ATACAGATAACAGATACATAGATAGACCGA";
	//std::string s1 = ref.substr(0, 25);
	//std::string s2 = ref.substr(5);
	std::string s1 = "AACAGAT";
	std::string s2 = "ACATGAA";

	StringVector input;
	input.push_back( s1 );
	input.push_back( s2 );


	std::cout << "S1: " << s1 << std::endl;
	std::cout << "S2: " << s2 << std::endl;
	//std::cout << "T:  " << t << std::endl;
	
	BWT b(s1);
	b.backwardSearch(s2);
	BWT c(s2);
	BWT d(input);
	d.printInfo();
	//std::string search = s1.substr(10);
	//std::cout << "Searching for " << search;
	//b.backwardSearch(search);
}

void stTest(std::string s)
{
	SuffixTree st(s + "$");
	st.printInfo();
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
			BWT b(line);
			b.printInfo();
			return;
		}

		if(count % 10000 == 0)
			std::cout << "Processed " << count << "\n";
		++count;
	}
	std::cout << "Loaded " << reads.size() << " reads\n";
	BWT b(reads);
}

//
//
//
int main(int /*argc*/, char** argv)
{
	//simpleTest();
	(void)argv;
	//dnaTest();
	//overlapTest();
	//stTest(argv[1]);
	//stReadTest();
	makeBWTFromReads(argv[1]);
	return 1;
}


