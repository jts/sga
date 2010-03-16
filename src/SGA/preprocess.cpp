//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// preprocess - prepare data files for assembly
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "preprocess.h"
#include "Timer.h"
#include "SeqReader.h"

//
// Getopt
//
#define SUBPROGRAM "preprocess"
static const char *PREPROCESS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *PREPROCESS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS1 READS2 ...\n"
"Prepare READS1, READS2, ... data files for assembly\n"
"Quality scores are assumed to be Sanger/Phred scaled - they must be converted if they are not\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the reads to FILE (default: basename(READS1).pp.fa)\n"
"      -q, --quality-trim=INT           perform Heng Li's BWT quality trim algorithm. \n"
"                                       Reads are trimmed according to the formula:\n"
"                                       argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT\n"
"                                       where l is the original read length.\n"
"      -m, --min-length=INT             discard sequences that are shorter than INT\n"
"                                       this is most useful when used in conjunction with --quality-trim\n"
"                                       note: this does not discard the pair of the read\n"
"      -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
"                                       the soft clip (quality-trim) option.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string outFile;
	static unsigned int qualityTrim = 0;
	static unsigned int hardClip = 0;
	static unsigned int minLength = 0;
	static bool bDiscardUncalled = true;
	static bool bIlluminaScaling = false;
}

static const char* shortopts = "o:q:m:h:vi";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "out",          required_argument, NULL, 'o' },
	{ "quality-trim", required_argument, NULL, 'q' },
	{ "hard-clip",    required_argument, NULL, 'h' },
	{ "min-length",   required_argument, NULL, 'm' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int preprocessMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("sga preprocess");
	parsePreprocessOptions(argc, argv);

	// If the output name hasn't been specified, make the outname
	// from the first file
	if(opt::outFile.empty())
	{
		assert(optind < argc);
		std::string firstReadsFile = argv[optind];
		std::string prefix = stripFilename(firstReadsFile);
		std::string suffix = getFileExtension(firstReadsFile);
		opt::outFile = prefix + ".pp." + suffix;
	}

	int64_t numReadsRead = 0;
	int64_t numReadsKept = 0;
	int64_t numBasesRead = 0;
	int64_t numBasesKept = 0;

	std::cerr << "Parameters:\n";
	std::cerr << "QualTrim: " << opt::qualityTrim << "\n";
	std::cerr << "HardClip: " << opt::hardClip << "\n";
	std::cerr << "Min length: " << opt::minLength << "\n";
	std::cerr << "Outfile: " << opt::outFile << "\n";
	if(opt::bDiscardUncalled)
		std::cerr << "Discarding sequences with uncalled bases\n";
	std::ofstream writer(opt::outFile.c_str());
	while(optind < argc)
	{
		std::string filename = argv[optind++];
		std::cerr << "Processing " << filename << "\n";
		SeqReader reader(filename);
		
		SeqRecord record;
		while(reader.get(record))
		{
			// Check if the sequence has uncalled bases
			std::string seqStr = record.seq.toString();
			std::string qualStr = record.qual;

			++numReadsRead;
			numBasesRead += seqStr.size();

			size_t nPos = seqStr.find_last_of('N');
			size_t dotPos = seqStr.find_last_of('.');

			bool bHasUncalled = nPos != std::string::npos || dotPos != std::string::npos;
			if(opt::bDiscardUncalled && bHasUncalled)
				continue; //skip read
			
			if(opt::qualityTrim > 0 && !qualStr.empty())
			{
				// Trim the sequence
				softClip(opt::qualityTrim, seqStr, qualStr);
				record.seq = seqStr;
				record.qual = qualStr;
			}
			else if(opt::hardClip > 0)
			{
				record.seq = seqStr.substr(0, opt::hardClip);
				if(!qualStr.empty())
					record.qual = qualStr.substr(0, opt::hardClip);
			}

			if(record.seq.length() < opt::minLength)
				continue;

			record.write(writer);

			++numReadsKept;
			numBasesKept += record.seq.length();
		}
	}
	
	writer.close();

	std::cerr << "Preprocess stats:\n";
	std::cerr << "Reads parsed:\t" << numReadsRead << "\n";
	std::cerr << "Reads kept:\t" << numReadsKept << " (" << (double)numReadsKept / (double)numReadsRead << ")\n"; 
	std::cerr << "Bases parsed:\t" << numBasesRead << "\n";
	std::cerr << "Bases kept:\t" << numBasesKept << " (" << (double)numBasesKept / (double)numBasesRead << ")\n"; 
	delete pTimer;
	return 0;
}

// Perform a soft-clipping of the sequence by removing low quality bases from the 
// 3' end using Heng Li's algorithm from bwa
void softClip(int qualTrim, std::string& seq, std::string& qual)
{
	assert(seq.size() == qual.size());

	int endpoint = 0; // not inclusive
	int max = 0;
	int i = seq.length() - 1;
	int terminalScore = Quality::char2phred(qual[i]);
	// Only perform soft-clipping if the last base has qual less than qualTrim
	if(terminalScore >= qualTrim)
		return;

	int subSum = 0;
	while(i >= 0)
	{
		int ps = Quality::char2phred(qual[i]);
		int score = qualTrim - ps;
		subSum += score;
		if(subSum > max)
		{
			max = subSum;
			endpoint = i;
		}
		--i;
	}

	// Clip the read
	seq = seq.substr(0, endpoint);
	qual = qual.substr(0, endpoint);
}

// 
// Handle command line arguments
//
void parsePreprocessOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'o': arg >> opt::outFile; break;
			case 'q': arg >> opt::qualityTrim; break;
			case 'i': arg >> opt::bIlluminaScaling; break;
			case 'm': arg >> opt::minLength; break;
			case 'h': arg >> opt::hardClip; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << PREPROCESS_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << PREPROCESS_VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 1) 
	{
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	} 

	if (die) 
	{
		std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
}
