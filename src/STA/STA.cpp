#include "STCommon.h"
#include "BWT.h"
#include "SuffixTree.h"

//
//
//
void simpleTest()
{
	std::string t = "GOOGOL$";
	std::cout << "T: " << t << "\n";
	BWT b(t);
	b.backwardSearch("GO");
}

//
//
//
void dnaTest()
{
	std::string ref = "ATACAGATAACAGATACATAGATAGACCGA";
	std::string s1 = ref.substr(0, 25);
	std::string s2 = ref.substr(5);
	std::string t = s1 + "$" + s2 + "&";

	std::cout << "S1: " << s1 << std::endl;
	std::cout << "S2:      " << s2 << std::endl;
	std::cout << "T:  " << t << std::endl;
	BWT b(t);
	std::string search = s1.substr(10);
	std::cout << "Searching for " << search;
	b.backwardSearch(search);
}

//
//
//
void overlapTest()
{
	std::string ref = "ATACAGATAACAGATACATAGATAGACCGA";
	SAStringVector input;
	input.push_back( SAString("str1", 0, ref.substr(0, 25) + "$") );
	input.push_back( SAString("str2", 0, ref.substr(5) + "&") );
	printVector(input);

	BWT b(input);
	SAString read = input.front();
	std::string s = read.str.substr(0, read.str.length() - 1);
	b.getOverlaps(s, 10);
}

void stTest(std::string s)
{
	SuffixTree st(s + "$");
	st.printInfo();
}

void stReadTest()
{
	std::string ref = "ATACAGATAACAGATACATAGATAGACCGAATAGACAGTACAGATACAGATATAGATAGATATAGACATAGA";
	SuffixTree sd(ref.substr(0, 25) + "$");

	for(int i = 0; i < 10; ++i)
	{
		sd.insert(ref.substr(i, 25) + "$");
		sd.printInfo();
	}
	/*
	std::string r1 = ref.substr(0,25) + "$";
	std::string r2 = ref.substr(1,25) + "$";
	std::string r3 = ref.substr(2,25) + "$";
	std::string r4 = ref.substr(3,25) + "$";

	SuffixTree sd(r1);
	sd.printInfo();
	sd.insert(r2);
	sd.printInfo();
	sd.insert(r3);
	sd.printInfo();
	sd.insert(r4);
	sd.printInfo();
	*/
}

//
//
//
int main(int /*argc*/, char** argv)
{
	//simpleTest();
	(void)argv;
	//overlapTest();
	//stTest(argv[1]);
	stReadTest();
	return 1;
}


