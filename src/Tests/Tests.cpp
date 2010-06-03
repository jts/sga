#include <iostream>
#include "Edge.h"
#include "Vertex.h"
#include "Bigraph.h"
#include "BWT.h"
#include "RLBWT.h"
#include "BWTWriter.h"

void dnaStringTests();

int main(int argc, char** argv)
{
    (void)argc;
    (void)argv;

    std::string file = argv[1];
    BWT* pBWT = new BWT(file);
    RLBWT* pRLBWT = new RLBWT(file);

    std::cout << "Standard BWT info:\n";
    pBWT->printInfo();

    std::cout << "\nRun-length BWT info: \n";
    pRLBWT->printInfo();


    std::cout << "\nTesting random access to RLBWT\n";
    for(size_t i = 0; i < pBWT->getBWLen(); ++i)
    {
        char b = pBWT->getChar(i);
        char r = pRLBWT->getChar(i);

       // printf("RLBWT[%zu] = %c, BWT[%zu] = %c\n", i, r, i, b);

        if(r != b)
        {
            printf("Test failed: RLBWT[%zu] expected %c, got %c\n", i, b, r);
            assert(false);
        }
    }

    pBWT->write("old.bwt");
    pRLBWT->write("new.bwt");

    delete pBWT;
    delete pRLBWT;

    return 0;
}

void dnaStringTests()
{
    std::string s1 = "GACACT";
    DNAString d1(s1);
    assert(s1 == d1.toString());
    std::cout << "S1: " << s1 << "\n";
    std::cout << "D1: " << d1.toString() << "\n";
    d1.reverse();
    std::cout << "Reversed: " << d1.toString() << "\n";

    d1 = "GGGACC";
    std::cout << "d1: " << d1.toString() << "\n";
    DNAString d2(d1);
    std::cout << "d2: " << d2.toString() << "\n"; 

    d1 = d2;
    std::cout << "d1: " << d1.toString() << "\n";
    d1 = d1;

    for(size_t i = 0; i < 2*d1.length(); ++i)
    {
        std::cout << i << " suffix: " << d1.getSuffix(i) << "\n";
        std::cout << i << " sufstr: " << d1.getSuffixString(i) << "\n";
    }
}

