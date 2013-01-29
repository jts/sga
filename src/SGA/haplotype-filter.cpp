//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// haplotype-filter - perform quality checks
// on haplotypes and their associated variants
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "haplotype-filter.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"
#include "HapgenUtil.h"
#include "Stats.h"
#include "HashMap.h"

// Functions
bool applyFilter(const std::string& kmer, const BWTIndexSet& indices);
std::vector<size_t> getPopulationCoverageCount(const std::string& kmer, const BWTIndexSet& indices);
double calculateLogProbabilityNullMultinomial(const std::vector<size_t>& sample_count);
double calculateLogProbabilitySegregating(const std::vector<size_t>& sample_count);
double LM(const std::vector<size_t>& sample_count);

//
// Getopt
//
#define SUBPROGRAM "haplotype-filter"
static const char *HAPLOTYPE_FILTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *HAPLOTYPE_FILTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... HAPLOTYPE_FILE VCF_FILE\n"
"Remove haplotypes and their associated variants from a data set.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reads=FILE                 load the FM-index of the reads in FILE\n"
"      -o, --out-prefix=STR             write the passed haplotypes and variants to STR.vcf and STR.fa\n" 
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static size_t k = 41;
    static std::string readsFile;
    static std::string outFile = "hapfilter.vcf";
    static std::string haplotypeFile;
    static std::string vcfFile;
}

static const char* shortopts = "o:r:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "reads",                 required_argument, NULL, 'r' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int haplotypeFilterMain(int argc, char** argv)
{
    parseHaplotypeFilterOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    std::cout << "Loading FM-index of " << opt::readsFile << "...";
    BWTIndexSet indices;
    std::string prefix = stripGzippedExtension(opt::readsFile);
    indices.pBWT = new BWT(prefix + BWT_EXT);
    indices.pPopIdx = new PopulationIndex(prefix + POPIDX_EXT);
    indices.pSSA = new SampledSuffixArray(prefix + SAI_EXT, SSA_FT_SAI);
    indices.pQualityTable = new QualityTable();
    std::cout << "done\n";

    std::ofstream outFile(opt::outFile.c_str());
    std::ifstream inFile(opt::vcfFile.c_str());

    std::string line;
    while(getline(inFile, line))
    {
        assert(line.size() > 0);
        if(line[0] == '#')
        {
            outFile << line << "\n";
            continue;
        }

        StringVector fields = split(line, '\t');
        std::string kmer = fields[2];
        
        // HACK
        // Shrink the VCF kmer
        kmer = kmer.substr(10, opt::k);
        
        double val = LM(getPopulationCoverageCount(kmer, indices));
        std::stringstream lmss;
        lmss << fields[7];
        lmss << ";LM=" << val;
        fields[7] = lmss.str();
        for(size_t i = 0; i < fields.size(); ++i)
            outFile << fields[i] << "\t";
        outFile << "\n";
    }
    
    /*
    // Load a hash of k-mer -> haplotype
    StringVector haplotypes;
    HashMap<std::string, size_t> kmer_to_haplotype;

    // Load haplotypes into hash
    std::cout << "Loading haplotypes\n";
    size_t hash_k = 61; // HACK
    SeqReader reader(opt::haplotypeFile);
    SeqRecord record;
    while(reader.get(record))
    {
        haplotypes.push_back(record.seq.toString());
        std::string& haplotype = haplotypes.back();

        size_t nk = haplotypes.size() - hash_k + 1;
        for(size_t i = 0; i < nk; ++i)
        {
            std::string kmer = haplotype.substr(i, hash_k);
            kmer_to_haplotype[kmer] = haplotypes.size() - 1;
        }
//        applyFilter(record.seq.toString(), indices);
    }
    */

    // Cleanup
    delete indices.pBWT;
    delete indices.pPopIdx;
    delete indices.pQualityTable;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

//
bool applyFilter(const std::string& kmer, const BWTIndexSet& indices)
{
    printf("Kmer --- %s\n", kmer.c_str());
    std::vector<size_t> sample_coverage = getPopulationCoverageCount(kmer, indices);
    std::copy(sample_coverage.begin(), sample_coverage.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << "\n";

    double lr = LM(sample_coverage);
    return lr > 0;
}

// Calculate the probability of the count vector under the null using a multinomial distribution
double calculateLogProbabilityNullMultinomial(const std::vector<size_t>& sample_count)
{
    size_t num_samples = sample_count.size();
    size_t total_reads = 0;
    for(size_t i = 0; i < num_samples; ++i)
        total_reads += sample_count[i];
    double p = (double)total_reads / num_samples;
    double log_p = log(p);
    double sum = 0.0f;

    printf("N: %zu O: %zu P: %lf LP: %lf\n", num_samples, total_reads, p, log_p);

    for(size_t i = 0; i < num_samples; ++i)
    {
        size_t o_i = sample_count[i];
        double c = Stats::logPoisson(o_i, p);
        sum += c;
    }
    return sum;
}

double calculateLogProbabilitySegregating(const std::vector<size_t>& sample_count)
{
    /*
    size_t num_samples = sample_count.size();
    size_t total_reads = 0;
    size_t n_0 = 0; // number of samples with no depth 
    for(size_t i = 0; i < num_samples; ++i)
    {
        total_reads += sample_count[i];
        n_0 += sample_count[i] == 0 ? 1 : 0;
    }

    size_t d = 4; // HACK HACK, estimate this
    size_t M = std::max((size_t)round(total_reads/d), num_samples - n_0);
    double q = (double)total_reads / M;
    double r = (double)M / N;

    printf("Estimate of M: %zu q: %lf\n", M, q);

    // Calculate logQ(M)
    double Q_M = 1 - (1 - exp(-q)) * r;
    printf("
    return 0.0f;
    */
    (void)sample_count;
    return 0.0f;
}

double LM(const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    std::map<size_t, size_t> samples_with_count;
    double o = 0;
    for(size_t i = 0; i < N; ++i)
    {
        o += sample_count[i];
        samples_with_count[sample_count[i]] += 1;
    }

    double n_0 = samples_with_count[0]; // number of samples with no depth
    double d = 4; // HACK HACK, estimate this
    //double M = std::max(round(o/d), N - n_0);
    double M = round(o/d);
    // Clamp M at the number of samples
    M = std::min((double)N, M);

    double p = o / N;
    double q = o / M;
    double r = M / N;
    double inv_r = N / M;

    double Q_M = 1 - (1 - exp(-q)) * r;
    double log_Q_M = log(Q_M);
    
    printf("o: %lf M: %lf\n", o, M);
    printf("n0: %lf q: %lf r: %lf log(Q(M)): %lf\n", n_0, q, r, log_Q_M);
    double t1 = n_0 * (p + log_Q_M);

    double t2 = 0.0;
    for(size_t k = 1; k < N; ++k)
    {
        double nk = samples_with_count[k];
        double tmp = nk * ((k - 1) * log(inv_r) - (q - p));
        t2 += tmp;
    }

    printf("T1: %lf T2: %lf LM: %lf\n", t1, t2, t1 + t2);
    return t1 + t2;
}
//
std::vector<size_t> getPopulationCoverageCount(const std::string& kmer, const BWTIndexSet& indices)
{
    //
    StringVector kmers;
    kmers.push_back(kmer);
    kmers.push_back(reverseComplement(kmer));
    SeqRecordVector records;
    HapgenUtil::extractHaplotypeReads(kmers, indices, opt::k, false, -1, 100000, &records, NULL);

    printf("Found %zu reads\n", records.size());

    // Count the number of reads per samples using their index in the BWT
    std::vector<size_t> sample_counts(indices.pPopIdx->getNumSamples());

    for(size_t i = 0; i < records.size(); ++i)
    {
        std::stringstream parser(records[i].id.substr(4));
        size_t index;
        parser >> index;
        size_t s_idx = indices.pPopIdx->getSampleIndex(index);
        sample_counts[s_idx]++;
    }
    return sample_counts;
}

// 
// Handle command line arguments
//
void parseHaplotypeFilterOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'r': arg >> opt::readsFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << HAPLOTYPE_FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << HAPLOTYPE_FILTER_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    opt::haplotypeFile = argv[optind++];
    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << HAPLOTYPE_FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
