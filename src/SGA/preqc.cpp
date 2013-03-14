//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// pre-QC - Perform pre-assembly quality checks on a set
//          of reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <queue>
#include "Util.h"
#include "preqc.h"
#include "Timer.h"
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SGACommon.h"
#include "HashMap.h"
#include "KmerDistribution.h"
#include "KmerOverlaps.h"
#include "BloomFilter.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filestream.h"

#if HAVE_OPENMP
#include <omp.h>
#endif

// Typedefs
typedef rapidjson::PrettyWriter<rapidjson::FileStream> JSONWriter;

// Functions

// Compute the distribution of unipath lengths from the k-de Bruijn graph.
// A branch in the graph with coverage c will be ignored when 
//    c / max_c < coverage_ratio_threshold 
//    where max_c is the highest-coverage branch
// 
void unipath_length_distribution(JSONWriter* pWriter,
                                 const BWTIndexSet& index_set,
                                 size_t k,
                                 double coverage_ratio_threshold, 
                                 size_t n_samples);

//
// Getopt
//
#define SUBPROGRAM "preqc"
static const char *PREQC_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *PREQC_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Perform pre-assembly quality checks\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads (default: 1)\n"
"          --reference=FILE             use the reference FILE to calculate GC plot\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string prefix;
    static std::string readsFile;
    static std::string referenceFile;
}

static const char* shortopts = "p:d:t:o:k:n:b:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "threads",            required_argument, NULL, 't' },
    { "prefix",             required_argument, NULL, 'p' },
    { "sample-rate",        required_argument, NULL, 'd' },
    { "kmer-size",          required_argument, NULL, 'k' },
    { "num-reads",          required_argument, NULL, 'n' },
    { "branch-cutoff",      required_argument, NULL, 'b' },
    { "reference",          required_argument, NULL, OPT_REFERENCE },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Returns the valid SUFFIX neighbors of the kmer
// in the de Bruijn graph. The neighbors are encoded
// as the single base they add. A coverage threshold
// is applied to filter out low-coverage extensions.
std::string get_valid_dbg_neighbors_ratio(const std::string& kmer,
                                          const BWTIndexSet& index_set,
                                          double coverage_ratio_threshold)
{
    std::string out;
    AlphaCount64 counts = 
        BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, 
                                                              index_set.pBWT,
                                                              ED_SENSE,
                                                              index_set.pCache);
    
    if(!counts.hasDNAChar())
        return out; // no extensions

    char max_b = counts.getMaxDNABase();
    size_t max_c = counts.get(max_b);

    for(size_t j = 0; j < 4; ++j)
    {
        char b = "ACGT"[j];
        size_t c = counts.get(b);
        if(c > 0 && (double)c / max_c >= coverage_ratio_threshold)
            out.push_back(b);
    }
    return out;
}

std::string get_valid_dbg_neighbors_coverage_and_ratio(const std::string& kmer,
                                                       const BWTIndexSet& index_set,
                                                       size_t min_coverage,
                                                       double min_ratio)
{
    std::string out;
    AlphaCount64 counts = 
        BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, 
                                                              index_set.pBWT,
                                                              ED_SENSE,
                                                              index_set.pCache);
    
    if(!counts.hasDNAChar())
        return out; // no extensions

    char max_b = counts.getMaxDNABase();
    size_t max_c = counts.get(max_b);

    for(size_t j = 0; j < 4; ++j)
    {
        char b = "ACGT"[j];
        size_t c = counts.get(b);
        if(c >= min_coverage && (double)c / max_c >= min_ratio)
            out.push_back(b);
    }
    return out;    
}

//
void generate_unipath_length_data(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    pWriter->String("UnipathLength");
    pWriter->StartArray();
    for(size_t k = 16; k < 96; k += 5)
        unipath_length_distribution(pWriter, index_set, k, 0.9, 1000);
    pWriter->EndArray();
}

//
void generate_kmer_coverage(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    pWriter->String("KmerDistribution");
    pWriter->StartObject();
    size_t n_samples = 10000;
    size_t k = 51;
    
    pWriter->String("k");
    pWriter->Int(k);
    
    KmerDistribution kmerDistribution;
    for(size_t i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        int n = s.size();
        int nk = n - k + 1;
        for(int j = 0; j < nk; ++j)
        {
            std::string kmer = s.substr(j, k);
            int count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
            kmerDistribution.add(count);
        }
    }

    pWriter->String("distribution");
    pWriter->StartArray();
    int max = kmerDistribution.getCutoffForProportion(0.95f);
    std::vector<int> count_vector = kmerDistribution.toCountVector(max);
    for(size_t i = 1; i < count_vector.size(); ++i)
    {
        pWriter->StartObject();
        pWriter->String("kmer-depth");
        pWriter->Int(i);
        pWriter->String("count");
        pWriter->Int(count_vector[i]);
        pWriter->EndObject();
    }

    pWriter->EndArray();
    pWriter->EndObject();
}

//
void generate_position_of_first_error(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    size_t n_samples = 100000;
    size_t k = 41;
    size_t starting_count = 10;
    size_t min_count = 3;

    std::vector<size_t> position_count;
    std::vector<size_t> error_count;
    for(size_t i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        size_t nk = s.size() - k + 1;
        size_t first_kmer_count = 
            BWTAlgorithms::countSequenceOccurrences(s.substr(0, k), index_set.pBWT);

        // Skip reads with a weak starting kmer
        if(first_kmer_count < starting_count)
            continue;

        for(size_t j = 1; j < nk; ++j)
        {
            size_t kmer_count =
                BWTAlgorithms::countSequenceOccurrences(s.substr(j, k), index_set.pBWT);

            if(j >= position_count.size())
            {
                position_count.resize(j+1);
                error_count.resize(j+1);
            }

            position_count[j] += 1;
            if(kmer_count < min_count)
            {
                error_count[j] += 1;
                break;
            }
        }
    }

    pWriter->String("FirstErrorPosition");
    pWriter->StartObject();

    pWriter->String("indices");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(i);
    pWriter->EndArray();

    pWriter->String("base_count");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(position_count[i]);
    pWriter->EndArray();

    pWriter->String("error_count");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(error_count[i]);
    pWriter->EndArray();
    pWriter->EndObject();
}

//
void generate_errors_per_base(JSONWriter* pWriter, const BWTIndexSet& index_set)
{

    int n_samples = 100000;
    size_t k = 25;

    double max_error_rate = 0.95;
    size_t min_overlap = 50;
    
    std::vector<size_t> position_count;
    std::vector<size_t> error_count;

    Timer timer("test", true);
#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        //KmerOverlaps::retrieveMatches(s, k, min_overlap, max_error_rate, 2, index_set);
        (void)k;
        KmerOverlaps::approximateMatch(s, min_overlap, max_error_rate, 2, 200, index_set);

        /*
        MultipleAlignment ma = 
            KmerOverlaps::buildMultipleAlignment(s, k, min_overlap, max_error_rate, 2, index_set);

        // Skip when there is insufficient depth to classify errors
        size_t ma_rows = ma.getNumRows();
        if(ma_rows <= 1)
            continue;

        size_t ma_cols = ma.getNumColumns();
        size_t position = 0;
        for(size_t j = 0; j < ma_cols; ++j)
        {
            char s_symbol = ma.getSymbol(0, j);

            // Skip gaps
            if(s_symbol == '-' || s_symbol == '\0')
                continue;
            
            SymbolCountVector scv = ma.getSymbolCountVector(j);
            int s_symbol_count = 0;
            char max_symbol = 0;
            int max_count = 0;

            for(size_t k = 0; k < scv.size(); ++k)
            {
                if(scv[k].symbol == s_symbol)
                    s_symbol_count = scv[k].count;
                if(scv[k].count > max_count)
                {
                    max_count = scv[k].count;
                    max_symbol = scv[k].symbol;
                }
            }

            //printf("P: %zu S: %c M: %c MC: %d\n", position, s_symbol, max_symbol, max_count);

            // Call an error at this position if the consensus symbol differs from the read
            //    and the support for the read symbol is less than 4 and the consensus symbol
            //    is strongly supported.
            bool is_error = s_symbol != max_symbol && s_symbol_count < 4 && max_count >= 3;

#if HAVE_OPENMP
            #pragma omp critical
#endif
            {
                if(position >= position_count.size())
                {
                    position_count.resize(position+1);
                    error_count.resize(position+1);
                }

                position_count[position]++;
                error_count[position] += is_error;
            }
            position += 1;
        }
        */
    }
    
    pWriter->String("ErrorsPerBase");
    pWriter->StartObject();
    
    pWriter->String("base_count");
    pWriter->StartArray();
    for(size_t i = 0; i < position_count.size(); ++i)
        pWriter->Int(position_count[i]);
    pWriter->EndArray();
    
    pWriter->String("error_count");
    pWriter->StartArray();
    for(size_t i = 0; i < position_count.size(); ++i)
        pWriter->Int(error_count[i]);
    pWriter->EndArray();

    pWriter->EndObject();
}

// Generate local graph complexity measure
void generate_local_graph_complexity(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 50000;
    size_t min_coverage_to_test = 5;
    size_t min_coverage_for_branch = 3;
    double min_coverage_ratio = 0.5f;

    pWriter->String("LocalGraphComplexity");
    pWriter->StartArray();
    for(size_t k = 16; k < 86; k += 5)
    {
        size_t num_branches = 0;
        size_t num_kmers = 0;

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;
            
            for(size_t j = 0; j < s.size() - k + 1; ++j)
            {
                std::string kmer = s.substr(j, k);
                size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set);
                if(count < min_coverage_to_test)
                    break;

                std::string extensions = 
                    get_valid_dbg_neighbors_coverage_and_ratio(kmer, 
                                                               index_set, 
                                                               min_coverage_for_branch, 
                                                               min_coverage_ratio);
#if HAVE_OPENMP
                #pragma omp critical
#endif
                {    
                    num_branches += extensions.size() > 1;
                    num_kmers += 1;
                }

            }
        }

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("num_kmers");
        pWriter->Int(num_kmers);
        pWriter->String("num_branches");
        pWriter->Int(num_branches);
        pWriter->EndObject();
    }
    pWriter->EndArray();
}

// Generate random walk length
void generate_random_walk_length(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 1000;
    size_t min_coverage = 5;
    double coverage_cutoff = 0.75f;
    size_t max_length = 30000;
    
    // Create a bloom filter to mark
    // visited kmers. We do not allow a new
    // walk to start at one of these kmers
    size_t bf_overcommit = 20;
    BloomFilter* bloom_filter = new BloomFilter;;
    bloom_filter->initialize(n_samples * max_length * bf_overcommit, 3);

    pWriter->String("RandomWalkLength");
    pWriter->StartArray();

    for(size_t k = 16; k < 86; k += 5)
    {
        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("walk_lengths");
        pWriter->StartArray();
#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            size_t walk_length = 0;
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;
            std::string kmer = s.substr(0, k);
            if(bloom_filter->test(kmer.c_str(), k) || BWTAlgorithms::countSequenceOccurrences(kmer, index_set) < min_coverage)
                continue;
            bloom_filter->add(kmer.c_str(), k);

            while(walk_length < max_length) 
            {
                std::string extensions = get_valid_dbg_neighbors_ratio(kmer, index_set, coverage_cutoff);
                if(!extensions.empty())
                {
                    kmer.erase(0, 1);
                    kmer.append(1, extensions[rand() % extensions.size()]);
                    walk_length += 1;
                    bloom_filter->add(kmer.c_str(), k);
                }
                else
                {
                    break;
                }
            }
#if HAVE_OPENMP
        #pragma omp critical
#endif
            pWriter->Int(walk_length);
        }
        pWriter->EndArray();
        pWriter->EndObject();
    }

    pWriter->EndArray();
    delete bloom_filter;
}

void generate_gc_distribution(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 100000;
    size_t k = 31;
    double gc_bin_size = 0.05;

    std::vector<double> read_gc_sum;
    size_t read_gc_n = 0;
    
    std::vector<double> ref_gc_sum;
    size_t ref_gc_n = 0;

    std::vector<size_t> coverage_vector;
    std::vector<double> gc_vector;

    read_gc_sum.resize(1.0f / gc_bin_size + 1);
    ref_gc_sum.resize(1.0f / gc_bin_size + 1);

    // Calculate the gc content of sampled reads
    for(int i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        if(s.size() < k)
            continue;

        size_t cov = BWTAlgorithms::countSequenceOccurrences(s.substr(0, k), index_set.pBWT);

        double gc = 0.f;
        double at = 0.f;
        for(size_t j = 0; j < s.size(); ++j)
        {
            if(s[j] == 'C' || s[j] == 'G')
                gc += 1;
            else
                at += 1;
        }
        
        double gc_f = gc / (gc + at);
        size_t bin_idx = gc_f / gc_bin_size;
        read_gc_sum[bin_idx] += gc_f;
        read_gc_n += 1;

        gc_vector.push_back(gc_f);
        coverage_vector.push_back(cov);
    }

    // If a reference file is provided, calculate the reference GC for comparison
    if(!opt::referenceFile.empty())
    {
        size_t window_size = 100;
        SeqReader reader(opt::referenceFile, SRF_NO_VALIDATION);
        SeqRecord record;
        while(reader.get(record))
        {
            DNAString& seq = record.seq;
            size_t l = seq.length();
            for(size_t i = 0; i < l; i += window_size)
            {
                double gc = 0.f;
                double at = 0.f;
                for(size_t j = i; j < i + window_size && j < l; ++j)
                {
                    char c = seq.get(j);
                    if(c == 'C' || c == 'G')
                        gc += 1;
                    else if(c == 'A' || c == 'T') // IUPAC symbols may be found
                        at += 1;
                }
                double gc_f = gc / (gc + at);
                size_t bin_idx = gc_f / gc_bin_size;
                ref_gc_sum[bin_idx] += gc_f;
                ref_gc_n += 1;
            }
        }
    }

    pJSONWriter->String("GCDistribution");
    pJSONWriter->StartObject();

    pJSONWriter->String("gc_bins");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < read_gc_sum.size(); ++i)
        pJSONWriter->Double(i * gc_bin_size + gc_bin_size / 2);
    pJSONWriter->EndArray();

    pJSONWriter->String("read_gc_prop");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < read_gc_sum.size(); ++i)
        pJSONWriter->Double(read_gc_sum[i] / read_gc_n);
    pJSONWriter->EndArray();

    if(!opt::referenceFile.empty())
    {
        pJSONWriter->String("ref_gc_prop");
        pJSONWriter->StartArray();
        for(size_t i = 0; i < read_gc_sum.size(); ++i)
            pJSONWriter->Double(ref_gc_sum[i] / ref_gc_n);
        pJSONWriter->EndArray();
    }

    pJSONWriter->String("gc_samples");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < gc_vector.size(); ++i)
        pJSONWriter->Double(gc_vector[i]);
    pJSONWriter->EndArray();

    pJSONWriter->String("cov_samples");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < gc_vector.size(); ++i)
        pJSONWriter->Double(coverage_vector[i]);
    pJSONWriter->EndArray();

    pJSONWriter->EndObject();
}

void generate_duplication_rate(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 10000;
    size_t k = 50;

    size_t total_pairs = index_set.pBWT->getNumStrings() / 2;
    size_t num_pairs_checked = 0;
    size_t num_duplicates = 0;
#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        // Choose a read pair
        int64_t source_pair_idx = rand() % total_pairs;
        std::string r1 = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2);
        std::string r2 = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2 + 1);

        // Get the interval for $k1/$k2 which corresponds to the 
        // lexicographic rank of reads starting with those kmers
        std::string k1 = "$" + r1.substr(0, k);
        std::string k2 = "$" + r2.substr(0, k);

        BWTInterval i1 = BWTAlgorithms::findInterval(index_set.pBWT, k1);
        BWTInterval i2 = BWTAlgorithms::findInterval(index_set.pBWT, k2);

        std::vector<int64_t> pair_ids;
        for(int64_t j = i1.lower; j <= i1.upper; ++j)
        {
            int64_t read_id = index_set.pSSA->lookupLexoRank(j);
            if(read_id % 2 == 1)
                continue;
            
            int64_t pair_id = read_id % 2 == 0 ? read_id / 2 : (read_id - 1) / 2;
            if(pair_id != source_pair_idx)
                pair_ids.push_back(pair_id);
        }

        for(int64_t j = i2.lower; j <= i2.upper; ++j)
        {
            int64_t read_id = index_set.pSSA->lookupLexoRank(j);
            if(read_id % 2 == 0)
                continue;
            int64_t pair_id = read_id % 2 == 0 ? read_id / 2 : (read_id - 1) / 2;
            if(pair_id != source_pair_idx)
                pair_ids.push_back(pair_id);
        }

        std::sort(pair_ids.begin(), pair_ids.end());
        std::vector<int64_t>::iterator iter = 
            std::adjacent_find(pair_ids.begin(), pair_ids.end());
                                           
        bool has_duplicate = iter != pair_ids.end();
#if HAVE_OPENMP
        #pragma omp critical
#endif
        {
            num_pairs_checked += 1;
            num_duplicates += has_duplicate;
        }
    }

    pJSONWriter->String("PCRDuplicates");
    pJSONWriter->StartObject();
    pJSONWriter->String("num_duplicates");
    pJSONWriter->Int(num_duplicates);
    pJSONWriter->String("num_pairs");
    pJSONWriter->Int(num_pairs_checked);
    pJSONWriter->EndObject();
}

// Write a stream of calculated fragments sizes to the JSON file
void generate_pe_fragment_sizes(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 100000;
    size_t k = 51;
    size_t MAX_INSERT = 1500;

    std::vector<size_t> fragment_sizes;

    size_t total_pairs = index_set.pBWT->getNumStrings() / 2;
#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        // Choose a read pair
        int64_t source_pair_idx = rand() % total_pairs;
        std::string start_kmer = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2).substr(0, k);
        std::string end_kmer = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2 + 1).substr(0, k);
        // We assume that the pairs are orientated F/R therefore we reverse-complement k_end
        end_kmer = reverseComplement(end_kmer);

        // Aggressively walk the de Bruijn graph starting from k_start until k_end is found or we give up
        size_t steps = 0;
        bool found = false;
        while(!found && steps < MAX_INSERT)
        {
            // A coverage ratio of 1.0 will force use to only use the highest-coverage branch
            // This may generate erroneous insert sizes in (rare?) cases but will give a good approximation
            // to the real distribution
            std::string extensions = get_valid_dbg_neighbors_ratio(start_kmer, index_set, 1.0f);
            if(extensions.empty())
                break;

            start_kmer.erase(0, 1);
            start_kmer.append(1, extensions[0]);
            found = start_kmer == end_kmer;
            steps += 1;
        }

        if(found)
        {
#if HAVE_OPENMP
            #pragma omp critical
#endif
            fragment_sizes.push_back(steps + k);
        }
    }

    pJSONWriter->String("FragmentSize");
    pJSONWriter->StartObject();
    pJSONWriter->String("sizes");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < fragment_sizes.size(); ++i)
        pJSONWriter->Int(fragment_sizes[i]);
    pJSONWriter->EndArray();
    pJSONWriter->EndObject();
}

// Generate a report of the quality of each base
void generate_quality_stats(JSONWriter* pJSONWriter, const std::string& filename)
{
    size_t max_reads = 10000000;
    double sample_rate = 0.05;
    SeqReader reader(filename, SRF_KEEP_CASE | SRF_NO_VALIDATION);
    SeqRecord record;

    size_t n_reads = 0;
    std::vector<size_t> bases_checked;
    std::vector<size_t> sum_quality;
    std::vector<size_t> num_q30;

    while(reader.get(record) && n_reads++ < max_reads)
    {
        if((double)rand() / RAND_MAX < sample_rate && record.qual.length() == record.seq.length())
        {
            size_t l = record.seq.length();
            if(l > bases_checked.size())
            {
                bases_checked.resize(l);
                sum_quality.resize(l);
                num_q30.resize(l);
            }

            for(size_t i = 0; i < l; ++i)
            {
                bases_checked[i]++;
                size_t q = record.getPhredScore(i);
                sum_quality[i] += q;
                num_q30[i] += (q >= 30);
            }
        }
    }

    pJSONWriter->String("QualityScores");
    pJSONWriter->StartObject();
    
    pJSONWriter->String("mean_quality");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < bases_checked.size(); ++i)
        pJSONWriter->Double((float)sum_quality[i] / bases_checked[i]);
    pJSONWriter->EndArray();

    pJSONWriter->String("fraction_q30");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < bases_checked.size(); ++i)
        pJSONWriter->Double((float)num_q30[i] / bases_checked[i]);
    pJSONWriter->EndArray();
    pJSONWriter->EndObject();
}

// Main
//
int preQCMain(int argc, char** argv)
{
    parsePreQCOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    fprintf(stderr, "Loading FM-index of %s\n", opt::readsFile.c_str());
    BWTIndexSet index_set;
    index_set.pBWT = new BWT(opt::prefix + BWT_EXT);
    index_set.pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
    index_set.pCache = new BWTIntervalCache(10, index_set.pBWT);
    
    rapidjson::FileStream f(stdout);
    JSONWriter writer(f);

    // Top-level document
    writer.StartObject();

    generate_errors_per_base(&writer, index_set);

    /*
    generate_gc_distribution(&writer, index_set);
    generate_quality_stats(&writer, opt::readsFile);
    generate_pe_fragment_sizes(&writer, index_set);
    generate_kmer_coverage(&writer, index_set);
    generate_position_of_first_error(&writer, index_set);
    //generate_errors_per_base(&writer, index_set);
    generate_unipath_length_data(&writer, index_set);
    generate_duplication_rate(&writer, index_set);
    generate_random_walk_length(&writer, index_set);
    generate_local_graph_complexity(&writer, index_set);
    */

    // End document
    writer.EndObject();

    delete index_set.pBWT;
    delete index_set.pSSA;
    delete index_set.pCache;
    delete pTimer;
    return 0;
}

//
void unipath_length_distribution(JSONWriter* pWriter,
                                 const BWTIndexSet& index_set, 
                                 size_t k,
                                 double coverage_ratio_threshold, 
                                 size_t n_samples)
{
    pWriter->StartObject();
    pWriter->String("k");
    pWriter->Int(k);
    pWriter->String("walk_lengths");
    pWriter->StartArray();

#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < (int)n_samples; ++i)
    {
        // Get a random read from the BWT
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);

        // Use the first-kmer of the read to seed the seach
        if(s.size() < k)
            continue;
        
        HashMap<std::string, bool> loop_check;
        std::string start_kmer = s.substr(0, k);
        std::string curr_kmer = start_kmer;

        size_t walk_length = 0;
        bool done = false;
        while(!done)
        {
            loop_check[curr_kmer] = true;
            std::string extensions = get_valid_dbg_neighbors_ratio(curr_kmer, index_set, coverage_ratio_threshold);
            if(extensions.size() == 1)
            {
                curr_kmer.erase(0, 1);
                curr_kmer.append(1, extensions[0]);

                if(loop_check.find(curr_kmer) == loop_check.end())
                    walk_length += 1;
                else
                    done = true;
            }
            else
            {
                done = true;
            }
        }
#if HAVE_OPENMP
        #pragma omp critical
#endif
        pWriter->Int(walk_length);
    }
    pWriter->EndArray();
    pWriter->EndObject();
}
// 
// Handle command line arguments
//
void parsePreQCOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_HELP:
                std::cout << PREQC_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PREQC_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << PREQC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripGzippedExtension(opt::readsFile);
}
