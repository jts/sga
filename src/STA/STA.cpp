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
}


//
//
//
int main(int /*argc*/, char** argv)
{
	//simpleTest();
	(void)argv;
	//overlapTest();
	stTest(argv[1]);
	return 1;
}


