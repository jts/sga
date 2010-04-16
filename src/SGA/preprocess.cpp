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
"Quality scores are assumed to be Sanger/Phred scaled - they must be converted if they are not.\n"
"If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -o, --out=FILE                   write the reads to FILE (default: basename(READS1).pp.fa)\n"
"      -p, --pe-mode=INT                0 - do not treat reads as paired\n"
"                                       1 - reads are paired with the first read in READS1 and the second\n"
"                                       read in READS2. The paired reads will be interleaved in the output file\n"
"      -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm. \n"
"                                       Reads are trimmed according to the formula:\n"
"                                       argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT\n"
"                                       where l is the original read length.\n"
"      -m, --min-length=INT             discard sequences that are shorter than INT\n"
"                                       this is most useful when used in conjunction with --quality-trim\n"
"      -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
"                                       the soft clip (quality-trim) option.\n"
"      -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string outFile;
	static unsigned int qualityTrim = 0;
	static unsigned int hardClip = 0;
	static unsigned int minLength = 0;
	static unsigned int peMode = 0;
	static double sampleFreq = 1.0f;
	static bool bDiscardUncalled = true;
	static bool bIlluminaScaling = false;
}

static const char* shortopts = "o:q:m:h:p:s:vi";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "out",          required_argument, NULL, 'o' },
	{ "quality-trim", required_argument, NULL, 'q' },
	{ "pe-mode",      required_argument, NULL, 'p' },
	{ "hard-clip",    required_argument, NULL, 'h' },
	{ "min-length",   required_argument, NULL, 'm' },
	{ "sample",       required_argument, NULL, 's' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static int64_t s_numReadsRead = 0;
static int64_t s_numReadsKept = 0;
static int64_t s_numBasesRead = 0;
static int64_t s_numBasesKept = 0;

//
// Main
//
int preprocessMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("sga preprocess");
	parsePreprocessOptions(argc, argv);

	std::cerr << "Parameters:\n";
	std::cerr << "QualTrim: " << opt::qualityTrim << "\n";
	std::cerr << "HardClip: " << opt::hardClip << "\n";
	std::cerr << "Min length: " << opt::minLength << "\n";
	std::cerr << "Sample freq: " << opt::sampleFreq << "\n";
	std::cerr << "PE Mode: " << opt::peMode << "\n";
	std::cerr << "Outfile: " << (opt::outFile.empty() ? "stdout" : opt::outFile) << "\n";
	if(opt::bDiscardUncalled)
		std::cerr << "Discarding sequences with uncalled bases\n";

	std::ostream* pWriter;
	if(opt::outFile.empty())
	{
		pWriter = &std::cout;
	}
	else
	{
		std::ofstream* pFile = new std::ofstream(opt::outFile.c_str());
		assertFileOpen(*pFile, opt::outFile);
		pWriter = pFile;
	}

	if(opt::peMode == 0)
	{
		// Treat files as SE data
		while(optind < argc)
		{
			std::string filename = argv[optind++];
			std::cerr << "Processing " << filename << "\n";
			SeqReader reader(filename);
			SeqRecord record;

			while(reader.get(record))
			{
				bool passed = processRead(record);
				if(passed && samplePass())
				{
					record.write(*pWriter);
					++s_numReadsKept;
					s_numBasesKept += record.seq.length();
				}
			}
		}
	}
	else
	{
		assert(opt::peMode == 1);
		int numFiles = argc - optind;
		if(numFiles % 2 == 1)
		{
			std::cerr << "Error: An even number of files must be given for pe-mode 1\n";
			exit(EXIT_FAILURE);
		}

		while(optind < argc)
		{
			std::string filename1 = argv[optind++];
			std::string filename2 = argv[optind++];
			
			SeqReader reader1(filename1);
			SeqReader reader2(filename2);

			std::cerr << "Processing pe files" << filename1 << ", " << filename2 << "\n";
			SeqRecord record1;
			SeqRecord record2;
			while(reader1.get(record1) && reader2.get(record2))
			{
				// Ensure the read names are sensible
				std::string expectedID2 = getPairID(record1.id);
				std::string expectedID1 = getPairID(record2.id);

				if(expectedID1 != record1.id || expectedID2 != record2.id)
				{
					std::cerr << "Warning: Pair IDs do not match (expected format /1,/2 or /A,/B)\n";
				}

				bool passed1 = processRead(record1);
				bool passed2 = processRead(record2);

				if(passed1 && passed2 && samplePass())
				{
					record1.write(*pWriter);
					record2.write(*pWriter);
					s_numReadsKept += 2;
					s_numBasesKept += record1.seq.length();
					s_numBasesKept += record2.seq.length();
				}
			}
		}
	}

	if(pWriter != &std::cout)
		delete pWriter;

	std::cerr << "Preprocess stats:\n";
	std::cerr << "Reads parsed:\t" << s_numReadsRead << "\n";
	std::cerr << "Reads kept:\t" << s_numReadsKept << " (" << (double)s_numReadsKept / (double)s_numReadsRead << ")\n"; 
	std::cerr << "Bases parsed:\t" << s_numBasesRead << "\n";
	std::cerr << "Bases kept:\t" << s_numBasesKept << " (" << (double)s_numBasesKept / (double)s_numBasesRead << ")\n"; 
	delete pTimer;
	return 0;
}

// Process a single read by quality trimming, filtering
// returns true if the read should be kept
bool processRead(SeqRecord& record)
{
	// Check if the sequence has uncalled bases
	std::string seqStr = record.seq.toString();
	std::string qualStr = record.qual;

	++s_numReadsRead;
	s_numBasesRead += seqStr.size();

	size_t nPos = seqStr.find_last_of('N');
	size_t dotPos = seqStr.find_last_of('.');

	bool bHasUncalled = nPos != std::string::npos || dotPos != std::string::npos;
	if(opt::bDiscardUncalled && bHasUncalled)
		return false; //skip read
	
	// Hard clip
	if(opt::hardClip > 0)
	{
		record.seq = seqStr.substr(0, opt::hardClip);
		if(!qualStr.empty())
			record.qual = qualStr.substr(0, opt::hardClip);
	}

	// Quality trim
	if(opt::qualityTrim > 0 && !qualStr.empty())
	{
		softClip(opt::qualityTrim, seqStr, qualStr);
		record.seq = seqStr;
		record.qual = qualStr;
	}

	if(record.seq.length() == 0 || record.seq.length() < opt::minLength)
		return false;
	return true;
}

// return true if the random value is lower than the acceptance value
bool samplePass()
{
	if(opt::sampleFreq >= 1.0f)
		return true; // no sampling
	
	double r = rand() / (RAND_MAX + 1.0f);
	return r < opt::sampleFreq;
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
			case 'p': arg >> opt::peMode; break;
			case 's': arg >> opt::sampleFreq; break;
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

	if(opt::peMode > 1)
	{
		std::cerr << SUBPROGRAM ": error pe-mode must be 0 or 1 (found: " << opt::peMode << ")\n";
		exit(EXIT_FAILURE);
	}
}
