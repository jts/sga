
#include <kmer-count.h>
#include <iostream>
#include <stack>
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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... input.bwt\n"
"Extract k-mers from a BWT-index.\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display program version\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -x, --kmer-threshold=N           Output only kmers whose multiplicity is >=N. (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";



namespace opt {
    static std::string bwtFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    static int kmerLength = 27;
    static int kmerThreshold = 1;
}

static const char* shortopts = "d:k:x:";
enum { OPT_HELP = 1, OPT_VERSION };
static const struct option longopts[] = {
    { "sample-rate",           required_argument, NULL, 'd' },
    { "kmer-size",             required_argument, NULL, 'k' },
    { "kmer-threshold",        required_argument, NULL, 'x' },
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
            case 'x': arg >> opt::kmerThreshold; break;
            case OPT_HELP:
                std::cout << KMERCOUNT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERCOUNT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    } else if (argc - optind > 1)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(opt::kmerThreshold <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::bwtFile = argv[optind++];
}



//
// BWT Traversal algorithm
//

// Stack structure used in the depth first search of kmers
struct stack_elt_t
{
    size_t strsz;
    char bp;
    BWTInterval range;
    stack_elt_t(size_t strsz,char bp):strsz(strsz),bp(bp){}
    stack_elt_t(size_t strsz,char bp,const BWTInterval& range):strsz(strsz),bp(bp),range(range){}
};

// extract all kmers of a bwt
void traverse_kmer(const BWT* pBWT, unsigned int k, unsigned int kmerThreshold)
{
    std::stack< stack_elt_t > stack;
    std::string str; // string storing the current path

    // Intitialize the search with root elements
    for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
    {
        stack_elt_t e(str.size(),ALPHABET[i]);
        BWTAlgorithms::initInterval(e.range,e.bp,pBWT);
        if (e.range.isValid()) stack.push(e);
    }

    // Perform the kmer search
    while(!stack.empty())
    {
        stack_elt_t top = stack.top();
        stack.pop();
        str.resize(top.strsz);
        str.push_back(top.bp);
        if (str.length()>=k) {
            // we reach a suffix of size k, print it
            if (top.range.size()>=kmerThreshold)
            {
                std::cout << reverse(str) << '\t' << top.range.size() << std::endl;
            }
        } else
        {
            // not yet a suffix of size k, push next candidates
            for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
            {
                stack_elt_t e(str.size(),ALPHABET[i],top.range);
                BWTAlgorithms::updateInterval(e.range,e.bp,pBWT);
                if (e.range.isValid()) stack.push(e);
            }
        }
    }
}




//
// Main
//

int kmerCountMain(int argc, char** argv)
{
    parseKmerCountOptions(argc,argv);
    BWT* pBWT = new BWT(opt::bwtFile,opt::sampleRate);
    traverse_kmer(pBWT,opt::kmerLength,opt::kmerThreshold);
    return 0;
}



