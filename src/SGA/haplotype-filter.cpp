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
double LM(double d, const std::vector<size_t>& sample_count);
double LMHaploidNaive(double d, const std::vector<size_t>& sample_count);
double LMDiploid(const std::vector<size_t>& sample_count);
std::vector<size_t> simulateCoverageNull(double d, size_t o, size_t N);
std::vector<size_t> simulateCoverageHaploid(double d, size_t o, size_t N);
void runSimulation();

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
    runSimulation();
    exit(0);

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
        
        printf("Kmer --- %s\n", kmer.c_str());
        std::vector<size_t> sample_coverage = getPopulationCoverageCount(kmer, indices);
        std::copy(sample_coverage.begin(), sample_coverage.end(), std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << "\n";

        double val = LMDiploid(sample_coverage);
        std::stringstream lmss;
        lmss << fields[7];
        lmss << ";LM=" << val << ";";
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
void runSimulation()
{
    size_t trials = 2000;
    size_t N = 1000;
    double d = 4;

    size_t min_o = 5;
    size_t max_o = round(N * d);
    (void)min_o;
    (void)max_o;

    double fpr = 0.5;

    size_t n_tp = 0;
    size_t n_fp = 0;
    size_t n_tn = 0;
    size_t n_fn = 0;

    for(size_t i = 0; i < trials; ++i)
    {
        double p = (double)rand() / RAND_MAX;
        size_t o = (rand() % (max_o - min_o)) + min_o;
        std::vector<size_t> counts;
        bool is_fp = p < fpr;
        if(is_fp)
            counts = simulateCoverageNull(d, o, N);
        else
            counts = simulateCoverageHaploid(d, o, N);

        if(opt::verbose > 0)
        {
            std::copy(counts.begin(), counts.end(), std::ostream_iterator<size_t>(std::cout, " "));
            std::cout << "\n";
        }

        double s = LM(d, counts);

        if(is_fp && s < 0)
            n_tn += 1;
        else if(is_fp && s > 0)
            n_fp += 1;
        else if(!is_fp && s < 0)
            n_fn += 1;
        else
            n_tp += 1;
        printf("%zu\t%d\t%zu\t%lf\n", i, is_fp, o, s);
    }
    fprintf(stderr, "tp: %zu fp: %zu tn: %zu fn: %zu\n", n_tp, n_fp, n_tn, n_fn);
}

//
bool applyFilter(const std::string& kmer, const BWTIndexSet& indices)
{
    printf("Kmer --- %s\n", kmer.c_str());
    std::vector<size_t> sample_coverage = getPopulationCoverageCount(kmer, indices);
    std::copy(sample_coverage.begin(), sample_coverage.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << "\n";

    double d = 4;
    double lr = LM(d, sample_coverage);
    return lr > 0;
}

double LM(double d, const std::vector<size_t>& sample_count)
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
    
    double t1 = n_0 * (p + log_Q_M);

    double t2 = 0.0;
    for(size_t k = 1; k < N; ++k)
    {
        double nk = samples_with_count[k];
        double tmp = nk * ((k - 1) * log(inv_r) - (q - p));
        t2 += tmp;
    }

    if(opt::verbose > 0)
    {
        printf("o: %lf M: %lf\n", o, M);
        printf("n0: %lf q: %lf r: %lf log(Q(M)): %lf\n", n_0, q, r, log_Q_M);
        printf("T1: %lf T2: %lf LM: %lf\n", t1, t2, t1 + t2);
    }

    return t1 + t2;
}

double LMHaploidNaive(double d, const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    double o = 0;
    size_t n0 = 0;
    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
    }

    //double M = std::max(round(o/d), N - n_0);
    double M = round(o/d);
    // Clamp M at the number of samples
    M = std::min((double)N, M);

    double p = o / N;
    double q = o / M;
    
    double sum = 0;
    for(size_t i = 0; i < sample_count.size(); ++i)
    {
        size_t oi = sample_count[i];
        if(oi == 0)
        {
            //sum += p - q;
            double h1_zero = (M/N) * exp(-q) + ((N - M) / N);
            double log_h0_zero = -p;
            double log_pzero = log(h1_zero) - log_h0_zero;
            sum += log_pzero;
        }
        else
        {
            //sum += (log(M/N) + Stats::logPoisson(oi, q) - Stats::logPoisson(oi, p));
            sum += log(M/N) + oi * log(q) - oi * log(p) + p - q;
        }
    }

    if(opt::verbose > 0)
    {
        printf("o: %lf M: %lf N: %zu p: %lf q: %lf n0: %zu\n", o, M, N, p, q, n0);
    }

    return sum;
}

double LMDiploid(const std::vector<size_t>& sample_count)
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
    // n_0 is an estimate of the number of ref-ref HOMs. Calculate 
    double ref_ref_proportion = n_0 / N;

    // Proportion of ref alleles in the population
    double hw_p = sqrt(ref_ref_proportion);

    // Proportion of non-ref alleles in the population
    double hw_q = 1.0f - hw_p;
    double d = 4; // HACK HACK, estimate this
    (void)d;

    // estimated number of non-ref alleles in the population
    double M = round(N * hw_q);
    double TWO_N = 2*N; // total alleles 
    assert(M < TWO_N);

    // Mean depth under the error model
    double p = o / N;

    // Mean depth of alt allele under segregating model
    double q = o / M;

    printf("LMDiploid params -- n0: %lf hw_p: %lf hw_q: %lf o: %lf p: %lf q: %lf\n", n_0, hw_p, hw_q, o, p, q);

    double sum = 0.0f;
    for(size_t i = 0; i < N; ++i)
    {
        size_t oi = sample_count[i];

        // Seg model
        double LOG_POIS = Stats::logPoisson(oi, q);
        double T1 = pow((M / TWO_N), 2) * pow(2, oi) * exp(-q);
        double T2 = ((TWO_N - M) / N) * (M / TWO_N);
        double S1 = log(T1 + T2);

        // Null model
        double LOG_NULL = Stats::logPoisson(oi, p);
        sum += (LOG_POIS + S1 - LOG_NULL);
    }
    printf("LMDiploid: %lf\n", sum);
    return sum;
}

std::vector<size_t> simulateCoverageNull(double /*d*/, size_t o, size_t N)
{
    // Distribute o reads over the N individuals
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        // Choose a random index and give it the count
        counts[rand() % (size_t)N]++;
    }
    return counts;
}

std::vector<size_t> simulateCoverageHaploid(double d, size_t o, size_t N)
{
    double M = round(o / d);
    if(opt::verbose > 0)
        printf("Simulate haploid %zu reads across %lf samples\n", o, M);
    assert(M <= N);

    // Select M indices at random to receive the counts
    std::vector<size_t> indices(N);
    for(size_t i = 0; i < N; ++i)
        indices[i] = i;

    std::random_shuffle(indices.begin(), indices.end());

    // Distribute o reads over M individuals
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        size_t idx = indices[rand() % (size_t)M];
        counts[idx]++;
    }

    return counts;
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
