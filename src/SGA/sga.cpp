//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// sga - Main assembler driver program
//
#include <string>
#include <iostream>
#include "index.h" 
#include "overlap.h"
#include "assemble.h"
#include "oview.h"
#include "rmdup.h"
#include "preprocess.h"
#include "merge.h"
#include "correct.h"
#include "subgraph.h"
#include "scaffold.h"

#define PROGRAM_BIN "sga"
#define AUTHOR "Jared Simpson"

void printUsage();

int main(int argc, char** argv)
{
    if(argc <= 1)
        printUsage();
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help")
        {
            printUsage();
            return 0;
        }

        if(command == "preprocess")
            preprocessMain(argc - 1, argv + 1);
        else if(command == "index")
            indexMain(argc - 1, argv + 1);
        else if(command == "merge")
            mergeMain(argc - 1, argv + 1);
        else if(command == "rmdup")
            rmdupMain(argc - 1, argv + 1);
        else if(command == "overlap")
            overlapMain(argc - 1, argv + 1);
        else if(command == "correct")
            correctMain(argc - 1, argv + 1);
        else if(command == "assemble")
            assembleMain(argc - 1, argv + 1);
        else if(command == "subgraph")
            subgraphMain(argc - 1, argv + 1);
        else if(command == "oview")
            oviewMain(argc - 1, argv + 1);
        else if(command == "scaffold")
            scaffoldMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
    }

    return 0;
}

//
void printUsage()
{
    std::cout << "Program: "  << PACKAGE_NAME << "\n";
    std::cout << "Version: " << PACKAGE_VERSION << "\n";
    std::cout << "Contact: " << AUTHOR << " [" << PACKAGE_BUGREPORT << "]\n";
    std::cout << "Usage: " << PROGRAM_BIN << " <command> [options]\n\n";
    std::cout << "Commands:\n";
    std::cout << "           preprocess   filter and quality-trim reads\n";
    std::cout << "           index        build the BWT and FM-index for a set of reads\n";
    std::cout << "           merge        merge multiple BWT/FM-index files into a single index\n";
    std::cout << "           correct      correct sequencing errors in a set of reads\n";
    std::cout << "           rmdup        remove duplicated or identical reads from the data set\n";
    std::cout << "           overlap      compute overlaps between reads\n";
    std::cout << "           assemble     generate contigs\n";
    std::cout << "           scaffold     generate ordered sets of contigs using distance estimates\n";
    std::cout << "           oview        view overlap alignments\n";
    std::cout << "           subgraph     extract a subgraph from a graph\n";
    std::cout << "\n\n";
}
