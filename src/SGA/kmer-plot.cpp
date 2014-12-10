//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by James Gurtowski (gurtowsk@cshl.edu)
// Released under the GPL
//-----------------------------------------------
//
// kmer-plot - Print kmer counts from one or more bwt's along a given input sequence
//

#include <kmer-count.h>
#include <iostream>
#include <stack>
#include <memory>
#include <getopt.h>
#include <BWT.h>
#include <BWTInterval.h>
#include <BWTAlgorithms.h>
#include "SeqReader.h"

//
// Getopt
//

#define SUBPROGRAM "kmer-plot"
static const char *KMERPLOT_VERSION_MESSAGE = SUBPROGRAM " Version " PACKAGE_VERSION "\n";
static const char *KMERPLOT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] in.fa test1.bwt [test2.bwt]\n"
"Generate a table of the k-mers in in.fa, and optionaly count the number of time they appears in testX.bwt.\n"
"Output on stdout kmers and their counts on forward and reverse strand\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display program version\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


namespace opt 
{
  static std::string inputSequenceFile;
  static std::vector<std::string> bwtFiles;
  static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
  static int kmerLength = 27;
  static int intervalCacheLength = 10;
}

static const char* shortopts = "d:k:x:";
enum { OPT_HELP = 1, OPT_VERSION };
static const struct option longopts[] = {
    { "sample-rate",           required_argument, NULL, 'd' },
    { "kmer-size",             required_argument, NULL, 'k' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


void parseKmerPlotOptions(int argc, char** argv) 
{

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case OPT_HELP:
                std::cout << KMERPLOT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERPLOT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::kmerLength <= 0 || opt::kmerLength % 2 == 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero and odd\n";
        std::cout << "\n" << KMERPLOT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(optind >= argc)
    {
      std::cerr << SUBPROGRAM ": missing input sequence file" << std::endl;
      std::cout << "\n" << KMERPLOT_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
    
    opt::inputSequenceFile = argv[optind++];
    
    for(;optind<argc;++optind)
    {
        opt::bwtFiles.push_back(argv[optind]);
    }

    if (opt::bwtFiles.size() < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        std::cout << "\n" << KMERPLOT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}


//
// Main
//

int kmerPlotMain(int argc, char** argv)
{

    // parse command line arguments
    parseKmerPlotOptions(argc,argv);

    // allocate BWT objects
    std::vector<BWTIndexSet> bwtIndicies;
    BWTIndexSet tmpIdx;

    for(std::vector<std::string>::iterator it = opt::bwtFiles.begin();
        it != opt::bwtFiles.end(); ++it)
    {
      std::cerr << "Loading " << *it << std::endl;

      tmpIdx.pBWT = new BWT(*it, opt::sampleRate);
      tmpIdx.pCache = new BWTIntervalCache(opt::intervalCacheLength, tmpIdx.pBWT);
      
      bwtIndicies.push_back(tmpIdx);
    }



    //Init sequence reader
    SeqReader reader(opt::inputSequenceFile, SRF_NO_VALIDATION);
    SeqRecord record;
    std::string kmer, kmer_rc;
    std::cerr << "Processing " << opt::inputSequenceFile << std::endl;


    size_t kmer_idx;
    std::vector<BWTIndexSet>::iterator indexset_it;

    //read the sequences from the file and split into kmers
    while(reader.get(record))
    {
      for(kmer_idx=0; kmer_idx < record.seq.length()-opt::kmerLength+1; ++kmer_idx)
      {
        kmer = record.seq.substr(kmer_idx, opt::kmerLength);
        kmer_rc = reverseComplement(kmer);

        std::cout << record.id << "\t" << kmer_idx << "\t" << kmer;
          
        //print out kmer count for all the bwts
        for(indexset_it = bwtIndicies.begin(); 
            indexset_it != bwtIndicies.end(); 
            ++indexset_it)
        {
          std::cout << '\t' << BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, *indexset_it);
          std::cout << '\t' << BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer_rc, *indexset_it);
        }
        
        std::cout << std::endl;

      }
    }

    // clean memory
    for(indexset_it = bwtIndicies.begin(); 
        indexset_it != bwtIndicies.end();
        ++indexset_it)
    {
      delete (*indexset_it).pBWT;
      delete (*indexset_it).pCache;
    }


    return 0;
}



