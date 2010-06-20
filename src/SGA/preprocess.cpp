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
#include "PrimerScreen.h"

static unsigned int DEFAULT_MIN_LENGTH = 40;
static int LOW_QUALITY_PHRED_SCORE = 3;

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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] --quality-scale=STR READS1 READS2 ...\n"
"Prepare READS1, READS2, ... data files for assembly\n"
"If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
//"      --quality-scale=STR              Specify the quality scaling to use. This parameter is mandatory and acceptable\n"
//"                                       strings are none, sanger, illumina1.3, illumina1.5. It is extremely important\n"
//"                                       to set this correctly.\n"
"      -o, --out=FILE                   write the reads to FILE (default: basename(READS1).pp.fa)\n"
"      -p, --pe-mode=INT                0 - do not treat reads as paired\n"
"                                       1 - reads are paired with the first read in READS1 and the second\n"
"                                       read in READS2. The paired reads will be interleaved in the output file\n"
"      -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm. \n"
"                                       Reads are trimmed according to the formula:\n"
"                                       argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT\n"
"                                       where l is the original read length.\n"
"      -f, --quality-filter=INT         discard the read if it contains more than INT low-quality bases.\n"
"                                       Bases with phred score <= 3 are considered low quality. Default: no filtering.\n"
"                                       The filtering is applied after trimming so bases removed are not counted.\n"
"      -m, --min-length=INT             discard sequences that are shorter than INT\n"
"                                       this is most useful when used in conjunction with --quality-trim\n"
"      -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
"                                       the soft clip (quality-trim) option.\n"
"      --permuteN                       If a basecall is N, randomly change it to one of ACGT instead of discarding the read.\n"
"                                       The quality value (if present) is not changed.\n"
"      -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum QualityScaling
{
    QS_UNDEFINED,
    QS_NONE,
    QS_SANGER,
    QS_ILLUMINA_1_3,
    QS_ILLUMINA_1_5
};

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static unsigned int qualityTrim = 0;
    static unsigned int hardClip = 0;
    static unsigned int minLength = DEFAULT_MIN_LENGTH;
    static int qualityFilter = -1;
    static unsigned int peMode = 0;
    static double sampleFreq = 1.0f;

    static bool bDiscardUncalled = true;
    static QualityScaling qualityScale = QS_SANGER;

    static bool bFilterGC = false;
    static double minGC = 0.0f;
    static double maxGC = 1.0;
    static bool bIlluminaScaling = false;
}

static const char* shortopts = "o:q:m:h:p:s:f:vi";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PERMUTE, OPT_QSCALE, OPT_MINGC, OPT_MAXGC };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "out",           required_argument, NULL, 'o' },
    { "quality-trim",  required_argument, NULL, 'q' },
    { "quality-filter",required_argument, NULL, 'f' },
    { "pe-mode",       required_argument, NULL, 'p' },
    { "hard-clip",     required_argument, NULL, 'h' },
    { "min-length",    required_argument, NULL, 'm' },
    { "sample",        required_argument, NULL, 's' },
    { "quality-scale", required_argument, NULL, OPT_QSCALE},
    { "min-gc",        required_argument, NULL, OPT_MINGC},
    { "max-gc",        required_argument, NULL, OPT_MAXGC},
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "permuteN",      no_argument,       NULL, OPT_PERMUTE },
    { NULL, 0, NULL, 0 }
};

static int64_t s_numReadsRead = 0;
static int64_t s_numReadsKept = 0;
static int64_t s_numBasesRead = 0;
static int64_t s_numBasesKept = 0;
static int64_t s_numReadsPrimer = 0;

//
// Main
//
int preprocessMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga preprocess");
    parsePreprocessOptions(argc, argv);

    std::cerr << "Parameters:\n";
    std::cerr << "QualTrim: " << opt::qualityTrim << "\n";

    if(opt::qualityFilter >= 0)
        std::cerr << "QualFilter: at most " << opt::qualityFilter << " low quality bases\n";
    else
        std::cerr << "QualFilter: no filtering\n"; 

    std::cerr << "HardClip: " << opt::hardClip << "\n";
    std::cerr << "Min length: " << opt::minLength << "\n";
    std::cerr << "Sample freq: " << opt::sampleFreq << "\n";
    std::cerr << "PE Mode: " << opt::peMode << "\n";
    std::cerr << "Quality scaling: " << opt::qualityScale << "\n";
    std::cerr << "MinGC: " << opt::minGC << "\n";
    std::cerr << "MaxGC: " << opt::maxGC << "\n";
    std::cerr << "Outfile: " << (opt::outFile.empty() ? "stdout" : opt::outFile) << "\n";
    if(opt::bDiscardUncalled)
        std::cerr << "Discarding sequences with uncalled bases\n";
    
    // Seed the RNG
    srand(time(NULL));

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
    std::cerr << "Reads failed primer screen:\t" << s_numReadsPrimer << " (" << (double)s_numReadsPrimer / (double)s_numReadsRead << ")\n";
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

    // Check for uncalled bases
    for(size_t i = 0; i < seqStr.size(); ++i)
    {
        if(seqStr[i] == 'N' || seqStr[i] == '.')
        {
            if(!opt::bDiscardUncalled)
                seqStr[i] = randomBase();
        }
        
        // Ensure base is upper case
        seqStr[i] = toupper(seqStr[i]);
    }

    // Ensure sequence is entirely ACGT
    size_t pos = seqStr.find_first_not_of("ACGT");
    if(pos != std::string::npos)
        return false;
    
    // Hard clip
    if(opt::hardClip > 0)
    {
        seqStr = seqStr.substr(0, opt::hardClip);
        if(!qualStr.empty())
            qualStr = qualStr.substr(0, opt::hardClip);
    }

    // Quality trim
    if(opt::qualityTrim > 0 && !qualStr.empty())
        softClip(opt::qualityTrim, seqStr, qualStr);

    // Quality filter
    if(opt::qualityFilter >= 0 && !qualStr.empty())
    {
        int numLowQuality = countLowQuality(seqStr, qualStr);
        if(numLowQuality > opt::qualityFilter)
            return false;
    }

    // Filter by GC content
    if(opt::bFilterGC)
    {
        double gc = calcGC(seqStr);
        if(gc < opt::minGC || gc > opt::maxGC)
            return false;
    }

    // Primer screen
    bool containsPrimer = PrimerScreen::containsPrimer(seqStr);
    if(containsPrimer)
    {
        ++s_numReadsPrimer;
        return false;
    }

    record.seq = seqStr;
    record.qual = qualStr;



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

// Count the number of low quality bases in the read
int countLowQuality(const std::string& seq, const std::string& qual)
{
    assert(seq.size() == qual.size());

    int sum = 0;
    for(size_t i = 0; i < seq.length(); ++i)
    {
        int ps = Quality::char2phred(qual[i]);
        if(ps <= LOW_QUALITY_PHRED_SCORE)
            ++sum;
    }
    return sum;
}

double calcGC(const std::string& seq)
{
    double num_gc = 0.0f;
    double num_total = 0.0f;
    for(size_t i = 0; i < seq.size(); ++i)
    {
        if(seq[i] == 'C' || seq[i] == 'G')
            ++num_gc;
        ++num_total;
    }
    return num_gc / num_total;
}

// 
char randomBase()
{
    int i = rand() % 4;
    switch(i)
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            assert(false);
    }
    return 'A';        
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
            case 'f': arg >> opt::qualityFilter; break;
            case 'i': arg >> opt::bIlluminaScaling; break;
            case 'm': arg >> opt::minLength; break;
            case 'h': arg >> opt::hardClip; break;
            case 'p': arg >> opt::peMode; break;
            case 's': arg >> opt::sampleFreq; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_PERMUTE: opt::bDiscardUncalled = false; break;
            case OPT_MINGC: arg >> opt::minGC; opt::bFilterGC = true; break;
            case OPT_MAXGC: arg >> opt::maxGC; opt::bFilterGC = true; break;
            case OPT_QSCALE:
                if(arg.str() == "none")
                {
                    opt::qualityScale = QS_NONE;
                }
                else if(arg.str() == "sanger")
                {
                    opt::qualityScale = QS_SANGER;
                }
                else if(arg.str() == "illumina1.3")
                {
                    opt::qualityScale = QS_ILLUMINA_1_3;
                    assert(false && "not implemented");
                }
                else if(arg.str() == "illumina1.5")
                {
                    opt::qualityScale = QS_ILLUMINA_1_5;
                    assert(false && "not implemented");
                }
                else
                    std::cout << "Unknown quality string value: " << arg.str() << "\n";
                    std::cout << "Expected one of: none, sanger, illumina1.3, illumina1.5\n";
                    exit(EXIT_FAILURE);
                break;
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

    if(opt::qualityScale == QS_UNDEFINED)
    {
        std::cerr << SUBPROGRAM ": required parameter --quality-scale not found, please specify this parameter\n";
        exit(EXIT_FAILURE);
    }

    if(opt::minLength < DEFAULT_MIN_LENGTH)
    {
        std::cerr << SUBPROGRAM ": WARNING - it is suggested that the min read length is " << DEFAULT_MIN_LENGTH << "\n";
        std::cerr << SUBPROGRAM ": Using very short reads may considerably impact the performance\n";
    }

    if(opt::qualityScale == QS_NONE && opt::qualityTrim > 0)
    {
        std::cerr << SUBPROGRAM ": If --quality-trim is specified, --quality-scale cannot be none\n";
        exit(EXIT_FAILURE);
    }
}
