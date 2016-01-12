//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson/James Gurtowski (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Kmer-count - Counting kmers in BWT or Sequence Files
//

#include <kmer-count.h>
#include <iostream>
#include <stack>
#include <memory>
#include <getopt.h>
#include <BWT.h>
#include <BWTInterval.h>
#include <BWTAlgorithms.h>

//
// Getopt
//

#define SUBPROGRAM "kmer-count"
static const char *KMERCOUNT_VERSION_MESSAGE = SUBPROGRAM " Version " PACKAGE_VERSION "\n";
static const char *KMERCOUNT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] src.{bwt,fa,fq} [test1.bwt] [test2.bwt]\n"
"Generate a table of the k-mers in src.{bwt,fa,fq}, and optionally count the number of time they appears in testX.bwt.\n"
"Output on stdout the canonical kmers and their counts on forward and reverse strand if input is .bwt\n"
"If src is a sequence file output forward and reverse counts for each kmer in the file\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display program version\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -c, --cache-length=N             Cache Length for bwt lookups (default: 10)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";



namespace opt {
    static std::string inputSequenceFile;
    static std::vector<std::string> bwtFiles;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    static unsigned int kmerLength = 27;
    static int intervalCacheLength = 10;
}

static const char* shortopts = "d:k:c:x:";
enum { OPT_HELP = 1, OPT_VERSION };
static const struct option longopts[] = {
    { "sample-rate",           required_argument, NULL, 'd' },
    { "kmer-size",             required_argument, NULL, 'k' },
    { "cache-length",          required_argument, NULL, 'c' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


void parseKmerCountOptions(int argc, char** argv) {

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'c': arg >> opt::intervalCacheLength; break;
            case OPT_HELP:
                std::cout << KMERCOUNT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERCOUNT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::kmerLength <= 0 || opt::kmerLength % 2 == 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero and odd\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(optind >= argc)
    {
      std::cerr << SUBPROGRAM ": missing input bwt/sequence file" << std::endl;
      std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }

    std::string inFile(argv[optind]);
    
    if(inFile.substr(inFile.find_last_of(".")) != ".bwt")
    {
      opt::inputSequenceFile = inFile;
      ++optind;
    }

    for(;optind<argc;++optind)
    {
        opt::bwtFiles.push_back(argv[optind]);
    }

    if (opt::inputSequenceFile.empty() && opt::bwtFiles.size() < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}



//
// BWT Traversal algorithm
//

// Stack structure used in the depth first search of kmers
struct stack_elt_t
{
    size_t str_sz;
    char bp;
    BWTInterval range;
    stack_elt_t(size_t strsz,char bp):str_sz(strsz),bp(bp){}
    stack_elt_t(size_t strsz,char bp,const BWTInterval& range):str_sz(strsz),bp(bp),range(range){}
};



// extract all canonical kmers of a bwt by performing a backward depth-first-search
void traverse_kmer(const std::vector<BWTIndexSet> &indicies)
{
    std::stack< stack_elt_t > stack;
    std::string str; // string storing the current path

    // Intitialize the search with root elements
    for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
    {
        stack_elt_t e(str.size(),ALPHABET[i]);
        BWTAlgorithms::initInterval(e.range,e.bp,indicies[0].pBWT);
        if (e.range.isValid()) stack.push(e);
    }

    std::vector<BWTIndexSet>::const_iterator idx_it;

    // Perform the kmer search
    while(!stack.empty())
    {
        // pop an element from the stack, and update the string path accordingly
        stack_elt_t top = stack.top();
        stack.pop();
        str.resize(top.str_sz);
        str.push_back(top.bp);
        if (str.length()>= opt::kmerLength) {
            // we found a kmer, retreive the number of occurence
            std::string seq(reverse(str));
            int64_t seq_count = top.range.size();

            // look for the count of its reverse complement
            std::string seq_rc(complement(str));
            BWTInterval range_rc = BWTAlgorithms::findInterval(indicies[0].pBWT,seq_rc);
            int64_t seq_rc_count = range_rc.isValid()?range_rc.size():0;

            // print the current kmer if canonical
            if (seq<seq_rc) {
              std::cout << seq << '\t' << seq_count << '\t' << seq_rc_count;
              
              for(idx_it = indicies.begin() + 1;
                  idx_it != indicies.end(); ++idx_it)
              {
                std::cout << "\t" << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seq, *idx_it);
                std::cout << "\t" << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seq_rc, *idx_it);
              }
              std::cout << std::endl;

            } else if (seq_rc_count<=0) {
                // the current kmer is not canonical, but the reverse complement doesn't exists
                // so print it now as it will never be traversed by the searching algorithm
                std::cout << seq_rc << '\t' << seq_rc_count << '\t' << seq_count;
                for(idx_it = indicies.begin() + 1;
                    idx_it != indicies.end(); ++idx_it)
                {
                  std::cout << "\t" << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seq_rc, *idx_it);
                  std::cout << "\t" << BWTAlgorithms::countSequenceOccurrencesSingleStrand(seq, *idx_it);
                }
                std::cout << std::endl;
            }
        } else
        {
            // not yet a suffix of size k, push next candidates
            for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
            {
                stack_elt_t e(str.size(),ALPHABET[i],top.range);
                BWTAlgorithms::updateInterval(e.range,e.bp,indicies[0].pBWT);
                if (e.range.isValid()) stack.push(e);
            }
        }
    }
}




void kmers_from_file(const std::string &inputFile,
                     const std::vector<BWTIndexSet> bwtIndicies)
{

    //Init sequence reader
    SeqReader reader(inputFile, SRF_NO_VALIDATION);
    SeqRecord record;
    std::string kmer, kmer_rc;
    std::cerr << "Processing " << inputFile << std::endl;

    size_t kmer_idx;
    std::vector<BWTIndexSet>::const_iterator indexset_it;

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
}

//
// Main
//

int kmerCountMain(int argc, char** argv)
{
    // parse command line arguments
    parseKmerCountOptions(argc,argv);

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

    if( opt::inputSequenceFile.empty() )
    {
      // run kmer search
      traverse_kmer(bwtIndicies);
    }else
    {
      kmers_from_file(opt::inputSequenceFile, bwtIndicies);
    }

    // clean memory
    std::vector<BWTIndexSet>::const_iterator indexset_it;
    for(indexset_it = bwtIndicies.begin(); 
        indexset_it != bwtIndicies.end();
        ++indexset_it)
    {
      delete (*indexset_it).pBWT;
      delete (*indexset_it).pCache;
    }
    return 0;
}



