//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// sga - Main assembler driver program
//
#include <string>
#include <iostream>
#include "index.h" 
#include "overlap.h"

#define VERSION "0.1"
#define PROGRAM_BIN "sga"
#define CONTACT "Jared Simpson [js18@sanger.ac.uk]"

void printUsage();

int main(int argc, char** argv)
{
	(void)argc;
	(void)argv;

	if(argc <= 1)
		printUsage();
	else
	{
		std::string command(argv[1]);
		if(command == "help" || command == "--help")
			printUsage();
		if(command == "index")
			indexMain(argc - 2, argv + 2);
		else if(command == "overlap")
			overlapMain(argc - 2, argv + 2);
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
	std::cout << "Program: "  << PROGRAM_BIN << "\n";
	std::cout << "Version: " << VERSION << "\n";
	std::cout << "Contact: " << CONTACT << "\n";
	std::cout << "Usage: " << PROGRAM_BIN << " <command> [options]\n\n";
	std::cout << "Command: index        index reads\n";
	std::cout << "         overlap      compute overlaps between reads\n";
	std::cout << "         assemble     generate contigs\n";
	std::cout << "\n\n";
}
