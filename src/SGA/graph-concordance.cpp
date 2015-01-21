//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// graph-concordance - check variants in a VCF file
// against an assembly graph
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <list>
#include "Util.h"
#include "graph-concordance.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"
#include "HapgenUtil.h"
#include "SGAStats.h"
#include "HashMap.h"
#include "DindelRealignWindow.h"
#include "VariantIndex.h"
#include "HapgenUtil.h"
#include "KmerDistribution.h"

//
// Getopt
//
#define SUBPROGRAM "graph-concordance"
static const char *GRAPH_CONCORDANCE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *GRAPH_CONCORDANCE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... -r variant.fastq -b base.fastq VCF_FILE\n"
"Count read support for variants in a vcf file\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"          --reference=STR              load the reference genome from FILE\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -g, --germline=FILE              load germline variants from FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int k = 51;
    static std::string variantFile;
    static std::string baseFile;
    static std::string haplotypeFile;
    static std::string vcfFile;
    static std::string referenceFile;
    static std::string germlineFile;
}

static const char* shortopts = "b:r:o:k:t:g:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_HAPLOID };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "reads",                 required_argument, NULL, 'r' },
    { "base",                  required_argument, NULL, 'b' },
    { "variant",               required_argument, NULL, 'r' },
    { "germline",              required_argument, NULL, 'g' },
    { "reference",             required_argument, NULL, OPT_REFERENCE },
    { "haploid",               no_argument,       NULL, OPT_HAPLOID },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Code ported from preqc for testing/development
//

// Enums
enum BranchClassification
{
    BC_NO_CALL,
    BC_ERROR,
    BC_VARIANT,
    BC_REPEAT
};

// struct to store some estimated properties of the genome
struct GenomeEstimates
{
    double average_read_length;
    size_t genome_size;
    size_t total_reads;

    // the mean number of reads starting at a given reference position
    double mean_read_starts;
};

// struct to store the parameters used by our classifier
struct ModelParameters
{
    int k;

    // the modal k-mer depth inferred from the count distribution
    double mode;

    // the kmer count model is a mixture of poissons representing
    // the copy number of a k-mer. this vector gives the mixture
    // proportions:
    // [ error_proportion, het_proportion, diploid_prop, 1-repeat, 2-repeat, ...]
    static const size_t NUM_MIXTURE_STATES = 10;
    static const size_t MAX_REPEAT_COPY = NUM_MIXTURE_STATES - 1;

    double mixture_proportions[NUM_MIXTURE_STATES];
    double mixture_means[NUM_MIXTURE_STATES];

    // This vector contains normalized proportions
    // for the repeat states. The indexing is the same
    // as the above arrays - for non-repeat states the
    // value will be 0.
    double repeat_mixture_normalized[NUM_MIXTURE_STATES];

    // prior on the per-base error rate
    static const double error_rate_prior;
};

// struct to store posterior probabilities calculated
// by the model
struct ModelPosteriors
{
    ModelPosteriors() : posterior_error(0.0f), 
                        posterior_variant(0.0f), 
                        posterior_repeat(0.0f),
                        classification(BC_NO_CALL) {}

    double posterior_error;
    double posterior_variant;
    double posterior_repeat;
    BranchClassification classification;
};

// Run the probablisitic classifier on a two-edge branch in the de Bruijn graph
extern ModelPosteriors classify_2_branch(const ModelParameters& params,
                                  const GenomeEstimates& genome,
                                  size_t higher_count, 
                                  size_t lower_count, 
                                  int delta,
                                  bool is_somatic);

// struct to store a description of the neighbors of a vertex in the graph
struct KmerNeighbors
{
    // This stores the occurrence count for each neighbor
    // including both strands.
    AlphaCount64 total_count;

    // Same as above but restricted to the kmers that have
    // coverage on both strands
    AlphaCount64 count_both_strands;

    // This stores the bases that have coverage on both strands
    std::string extensions_both_strands;

    // Return an ordering of the 4 possible extensions by their count
    // For instance A:2 C:5 G:8 T:0 will return "GCAT"
    static std::string getExtensionsFromCount(const AlphaCount64& a)
    {
        char sorted_bases[5] = "ACGT";
        sorted_bases[4] = '\0';
        a.getSorted(sorted_bases, 5);
        std::string r(sorted_bases);
        return r;
    }
};

extern int calculateHomopolymerLength(const VCFRecord& record, const ReadTable* refTable);


extern void learn_mixture_parameters(const KmerDistribution& distribution, ModelParameters& params);
extern size_t find_single_copy_peak(const KmerDistribution& distribution);

extern KmerNeighbors calculate_neighbor_data(const std::string& kmer, const BWTIndexSet& index_set);
extern KmerDistribution sample_kmer_counts(size_t k, size_t n, const BWTIndexSet& index_set);
extern int calculate_delta(const std::string& kmer, const KmerNeighbors nd, const BWTIndexSet& index_set);

extern ModelParameters calculate_model_parameters(size_t k, KmerDistribution& distribution);
extern ModelParameters calculate_model_parameters(size_t k, size_t samples, const BWTIndexSet& index_set);
extern GenomeEstimates estimate_genome_size_from_k_counts(size_t k, const BWTIndexSet& index_set);

bool operator==(const BWTIndexSet&bw1, const BWTIndexSet&bw2) {
    return bw1.pBWT == bw2.pBWT;
}

struct LRU {
    typedef std::pair<std::string,BWTIndexSet> keytype;

    struct keyCompare
    {
          bool operator() (const keytype& lhs, const keytype& rhs)
          {
              // sort by string first, using BWT index pointers to break ties
              if (lhs.first == rhs.first)
                  return lhs.second.pBWT < rhs.second.pBWT;
              else
                  return lhs.first < rhs.first;
          }
    };

    typedef std::list< keytype > timestamplisttype;
    typedef size_t valtype;
    typedef std::pair< valtype, timestamplisttype::iterator > valpairtype;

    static const size_t maxcount = 1000;

    typedef std::map< keytype, valpairtype, keyCompare> maptype ;

    maptype cache;
    timestamplisttype timestamplist;

    maptype::iterator locate( const std::string &s, const BWTIndexSet& index_set ) {
        keytype key( s, index_set );
        return cache.find( key );
    }

    void updateTimestamp( maptype::iterator& it ) {
        timestamplisttype::iterator lit = (it->second).second;
        keytype key = *lit;
        timestamplist.erase(lit);
        timestamplist.push_front(key);
        it->second.second = timestamplist.begin();
    }

    int evictlast() {
        if (timestamplist.size() == 0)
            return 1;
        keytype key = timestamplist.back();
        timestamplist.pop_back();
        maptype::iterator it = cache.find( key );
        cache.erase( it );
        return 0;
    }

    int findval( const std::string &s, const BWTIndexSet& index_set, valtype& value ) {
        std::map< keytype, valpairtype>::iterator it = locate( s, index_set );
        if ( it == cache.end() )
            return 1;
        else {
            updateTimestamp(it);
            value = (it->second).first;
            return 0;
        }
    }

    void add( const std::string &s, const BWTIndexSet& index_set, valtype& value ) {
        maptype::iterator it = locate( s, index_set );
        if ( it != cache.end() ) {
            updateTimestamp( it );
            (it->second).first = value;
        } else {
            if ( timestamplist.size() == maxcount ) 
                evictlast();
            keytype key(s, index_set);
            timestamplist.push_front( key );
            valpairtype val( value, timestamplist.begin() );
            cache[key] = val;
        }
    }

    void validateStructure() {
        if ( cache.size() != timestamplist.size() )
            std::cerr << "map size = " << cache.size() << " list size = " << timestamplist.size() << std::endl;

        for (maptype::iterator it = cache.begin(); it != cache.end(); ++it ) {
            keytype key = it->first;
            keytype foundkey = *((it->second).second);
            if (key != foundkey) {
                std::cerr << "map/list inconsistency: string = " << key.first << std::endl;
            }
        }
    }

};

size_t kmerCount(const std::string& x, const BWTIndexSet& index_set )
{
    static LRU memos;
    static size_t nlookups=0;

    if ( x.empty() )
        return 0;
    else {
        nlookups++;
        if (nlookups % 1000 == 0)
            std::cerr << " kmerCount - nlookups = " << nlookups << std::endl;

        size_t count;
        if ( !memos.findval(x, index_set, count) ) {
            size_t count2 = BWTAlgorithms::countSequenceOccurrences(x, index_set);
            if (count != count2) {
                std::cerr << "Inconsistency in cache - input " << x << std::endl;
                std::cerr << "Got " << count << " real answer " << count2 << std::endl;
                count = count2;
            }
            memos.validateStructure();
            return count;
        } else {
            count = BWTAlgorithms::countSequenceOccurrences(x, index_set);
            memos.add(x, index_set, count);
            memos.validateStructure();
            return count;
        }
    }
}

size_t kmerCountWithOneOff(const std::string& x, const BWTIndexSet& index_set )
{
    if(x.empty())
        return 0;

    size_t direct = BWTAlgorithms::countSequenceOccurrences(x, index_set);

    // One-offs, except for the terminal bases
    size_t oneoff = 0;
    std::string t = x;

    for(size_t i = 1; i < x.size() - 1; ++i)
    {
        char tmp = t[i];
        for(int bi = 0; bi < DNA_ALPHABET::size; ++bi)
        {
            char b = DNA_ALPHABET::getBase(bi);
            if(b != tmp)
            {
                t[i] = b;
                size_t count = kmerCount(x, index_set);
                oneoff += count;
            }
        }
        t[i] = tmp;
    }
    return direct + oneoff;
}


std::string applyVariant(const std::string& in, int pos,
                         const std::string& ref, const std::string& alt,
                         std::vector<int>& coordinate_map, 
                         bool update_map, int variant_tag = -1)
{
    assert(variant_tag < 0);
    std::string out = in;

    // Search the coordinate map for the index of the ref string to be substituted
    size_t idx = 0;
    while(idx < coordinate_map.size())
    {
        if(coordinate_map[idx] == pos)
            break;
        idx += 1;
    }

    // If we could not find the position in the coordinate map
    // then the input variant is for a reference base that has been removed
    // already. Skip it.
    if(idx == coordinate_map.size())
        return "";

    /*
    std::cout << "In: " << in << "-\n";
    std::cout << "Current map: \n";
    for(size_t i = 0; i < coordinate_map.size(); ++i)
        std::cout << i << " " << coordinate_map[i] << "\n";
    std::cout << "\n";
    std::cout << "Replace idx: " << idx << "\n";
    std::cout << "Replace coordinate: " << coordinate_map[idx] << "\n";
    std::cout << "ref: " << ref << "\n";
    std::cout << "Alt: " << alt << "\n";
    */

    // Ensure that the reference string at the variant matches the expected
    if(out.substr(idx, ref.length()) != ref)
    {
        std::cerr << "Warning: REF does not match reference\n";
        return "";
    }

    // Replace the reference sequence with the alternative
    out.replace(idx, ref.length(), alt);

    if(update_map)
    {
        // shift coordinate map to account for indels

        // Make an iterator to the first position that is replaced
        std::vector<int>::iterator fi = coordinate_map.begin() + idx;

        // Iterator which is one past the last position to remove
        std::vector<int>::iterator li = fi + ref.size();
        
        // erase returns an iterator pointing to the position
        // following the last position removed, this is where
        // we want to insert 
        std::vector<int>::iterator ii = coordinate_map.erase(fi, li);

        // Insert new elements for the incoming bases
        // These elements can be tagged with an integer
        // by the calling function to indicate the sample
        // the new bases came from
        coordinate_map.insert(ii, alt.size(), variant_tag);
        
        /*
        std::cout << "Updated map: \n";
        for(size_t i = 0; i < coordinate_map.size(); ++i)
            std::cout << coordinate_map[i] << " ";
        std::cout << "\n";
        */
        assert(coordinate_map.size() == out.size());
    }
    return out;
}

struct ReconstructedHaplotype
{
    std::string haplotype;

    std::string left_base_kmer;
    std::string left_somatic_kmer;
    
    std::string right_base_kmer;
    std::string right_somatic_kmer;

    size_t unique_kmers;
};

// Extract somatic k-mers from the haplotype accordining to the tagged
// entries in the coordinate map
StringVector extractSomaticKmers(std::string base_haplotype, size_t k, 
                                 ReconstructedHaplotype& rc)
                                 
{
    StringVector kmers;
    rc.unique_kmers = 0;
    bool previous_kmer_somatic = false;
    for(size_t i = 0; i < rc.haplotype.size() - k + 1; ++i)
    {
        std::string kmer = rc.haplotype.substr(i, k);
        bool is_somatic = base_haplotype.find(kmer) == std::string::npos;
        
        if(is_somatic)
        {
            kmers.push_back(kmer);
            rc.unique_kmers += 1;
            if(kmers.size() == 1 && i > 0)
            {
                // This is the leftmost somatic kmer
                rc.left_base_kmer = rc.haplotype.substr(i - 1, k);
                rc.left_somatic_kmer = kmer;
            }
        }

        if(!is_somatic && previous_kmer_somatic)
        {
            assert(i > 0);
            // This is the rightmost somatic kmer (so far)
            rc.right_base_kmer = kmer;
            rc.right_somatic_kmer = rc.haplotype.substr(i - 1, k);
        }

        previous_kmer_somatic = is_somatic;
    }

    return kmers;
}

int reconstructHaplotype(const VCFRecord& somatic_record, 
                         const VariantRecordVector nearby_germline, 
                         const ReadTable& refTable,
                         const BWTIndexSet& variantIndex,
                         const BWTIndexSet& baseIndex,
                         size_t flanking_size,
                         ReconstructedHaplotype& out)
{
    // Extract the reference haplotype
    int eventLength = somatic_record.varStr.length();
    int zeroBasedPos = somatic_record.refPosition - 1;
    int start = zeroBasedPos - flanking_size - 1;
    if (start < 0) 
        return -1;
    int end = zeroBasedPos + eventLength + 2 * flanking_size;
    const SeqItem& chr = refTable.getRead(somatic_record.refName);
    if(end > (int)chr.seq.length())
        return -1;

    std::string reference_haplotype = chr.seq.substr(start, end - start);
    
    // Set up a vector mapping positions on the current haplotype
    // to original reference coordinates. This allows use to apply
    // the germline variants to the haplotype when the coordinates
    // shift due to indels
    std::vector<int> coordinate_map(reference_haplotype.size());
    for(int i = 0; i < (int)reference_haplotype.size(); ++i)
        coordinate_map[i] = i + start;

    // Do not attempt haplotype reconstruction if there are too many variants
    if(nearby_germline.size() > 4)
        return -1;
    
    int SOMATIC_TAG = -2;
    int GERMLINE_TAG = -3;

    // Iterate over all possible combinations of germline variants
    size_t total_haplotypes = 1 << nearby_germline.size(); // 2 ^ n

    size_t best_haplotype_kmers = 0;
    std::string best_haplotype = "";

    for(size_t haplotype_id = 0; haplotype_id < total_haplotypes; haplotype_id++)
    {
        std::string haplotype = reference_haplotype;
        std::vector<int> haplotype_coordinates = coordinate_map;

        // Iterate over the germline variants and check if the i-th bit
        // is set in the haplotype id. If so, apply the variant to the haplotype
        for(size_t var_id = 0; var_id < nearby_germline.size(); ++var_id)
        {
            if( (haplotype_id >> var_id) & 1)
            {
                std::string new_haplotype = 
                    applyVariant(haplotype, 
                                 nearby_germline[var_id].position - 1,
                                 nearby_germline[var_id].ref_sequence, 
                                 nearby_germline[var_id].alt_sequence,
                                 haplotype_coordinates, 
                                 true, 
                                 GERMLINE_TAG);
                haplotype = new_haplotype;

                if(haplotype.empty())
                    break;
            }
        }

        if(haplotype.empty())
            continue;

        ReconstructedHaplotype current;

        // Apply the somatic variant to this haplotype
        current.haplotype = 
            applyVariant(haplotype, 
                         somatic_record.refPosition - 1,
                         somatic_record.refStr, 
                         somatic_record.varStr,
                         haplotype_coordinates, 
                         true, 
                         SOMATIC_TAG);
    
        if(current.haplotype.empty())
            continue;
        
        StringVector kmers = extractSomaticKmers(haplotype, flanking_size, current);
        
        // Count the number of somatic kmers that are present in the variant reads
        for(size_t ki = 0; ki < kmers.size(); ++ki)
        {
            int cv = kmerCount(kmers[ki], variantIndex);
            int cb = kmerCount(kmers[ki], baseIndex);
            int ov = HapgenUtil::getMaximumOneEdit(kmers[ki], variantIndex);
            int ob = HapgenUtil::getMaximumOneEdit(kmers[ki], baseIndex);

            fprintf(stderr, "Hap[%zu][%zu], %s %d %d %d %d\n", haplotype_id, ki, kmers[ki].c_str(), cv, cb, ov, ob);
            //if(BWTAlgorithms::countSequenceOccurrences(kmers[ki], variantIndex) > 0)
            //    read_kmers += 1;
        }

        size_t read_kmers = BWTAlgorithms::countSequenceOccurrences(current.left_somatic_kmer, variantIndex) +
                            BWTAlgorithms::countSequenceOccurrences(current.right_somatic_kmer, variantIndex);

        if(read_kmers > best_haplotype_kmers || best_haplotype_kmers == 0)
        {
            out = current;
            best_haplotype_kmers = read_kmers;
        }
    }

    return 0;
}

std::string get_neighbor(const std::string& x, EdgeDir dir, BWTIndexSet indexSet)
{
    int best_count = 0;
    std::string best = "";
    for(int i = 0; i < 4; i++)
    {
        std::string tmp = x;

        if(dir == ED_SENSE)
            tmp[opt::k - 1] = ALPHABET[i];
        else
            tmp[0] = ALPHABET[i];

        if(tmp == x)
            continue;
        
        int count = BWTAlgorithms::countSequenceOccurrences(tmp, indexSet);
        if(count >= best_count)
        {
            best_count = count;
            best = tmp;
        }
    }
    return best;
}

int pickLocalKAndReconstruct(const VCFRecord& somatic_record, 
                         const VariantRecordVector nearby_germline, 
                         const ReadTable& refTable,
                         const BWTIndexSet& variantIndex,
                         const BWTIndexSet& baseIndex,
                         const size_t flanking_size,
                         const ModelParameters& variant_params,
                         const ModelParameters& base_params,
                         ReconstructedHaplotype& out, 
                         size_t& local_k) {
    
    local_k = flanking_size;
    size_t last_k  = flanking_size;

    int sum_y, sum_z;
    double coverage;

    int err;
    if ((err = reconstructHaplotype(somatic_record, nearby_germline, refTable, variantIndex, baseIndex, local_k, out)))
        return err;
    do 
    {
        std::string left_y = out.left_somatic_kmer;
        std::string left_z = get_neighbor(left_y, ED_SENSE, variantIndex);

        std::string right_y = out.right_somatic_kmer;
        std::string right_z = get_neighbor(right_y, ED_ANTISENSE, variantIndex);

        sum_y = kmerCount(left_y, variantIndex) + kmerCount(right_y, variantIndex);
        sum_z = kmerCount(left_z, variantIndex) + kmerCount(right_z, variantIndex);

        coverage = (base_params.mixture_means[1] + variant_params.mixture_means[1])*base_params.k/local_k;
        fprintf(stderr,"pick local k: k = %d, sum_y = %d, sum_z = %d, coverage=%lf\n", (int)local_k, sum_y, sum_z, coverage);
        if (sum_y + sum_z > 12*coverage) /* 4 kmers + 2sd each */
        {
            last_k = local_k;
            local_k += 2;
            err = reconstructHaplotype(somatic_record, nearby_germline, refTable, variantIndex, baseIndex, local_k, out);
        } 
    } while (!err && ( (sum_y + sum_z) > 12*coverage ) );

    // if we errored out, go back and re-run the last k; otherwise, we've
    // gotten the right reconstruction, and just return 0.
    if (err) 
    {
        local_k = last_k;
        fprintf(stderr,"pick local k: chose k = %d\n", (int)local_k);
        return reconstructHaplotype(somatic_record, nearby_germline, refTable, variantIndex, baseIndex, local_k, out);
    } 
    else 
    {
        return 0;
    }
}

int new_delta_calculation(const std::string& x,
                          const std::string& y,
                          const std::string& z,
                          BWTIndexSet index,
                          EdgeDir dir)
{
    int delta = BWTAlgorithms::countSequenceOccurrences(y, index) +
                BWTAlgorithms::countSequenceOccurrences(z, index);
    
    std::string lmer = x;
    if(dir == ED_SENSE)
    {
        std::string xy = x + y[y.size() - 1];
        std::string xz = x + z[z.size() - 1];
        delta -= (BWTAlgorithms::countSequenceOccurrences(xy, index) +
                  BWTAlgorithms::countSequenceOccurrences(xz, index));

    } 
    else
    {
        std::string yx = y[0] + x;
        std::string zx = z[0] + x;
        delta -= (BWTAlgorithms::countSequenceOccurrences(yx, index) +
                  BWTAlgorithms::countSequenceOccurrences(zx, index));
    }
    assert(delta >= 0);
    return delta;
}

double apply_branch_model(const std::string& x,
                          const std::string& y,
                          const std::string& z,
                          const BWTIndexSet& index,
                          const ModelParameters& params,
                          const GenomeEstimates& genome,
                          bool is_somatic)
{
    int c_x = kmerCount(x, index);
    int c_y = kmerCount(y, index);
    int c_z = kmerCount(z, index);
    fprintf(stderr, "Counts: x: %d y: %d z: %d\n", c_x, c_y, c_z);
    
    int delta = 0; //calculate_delta(x, v_neighbors, variantIndex);
    int c_min = std::min(c_y, c_z);
    int c_max = std::max(c_y, c_z);

    ModelPosteriors ret = classify_2_branch(params, genome, c_max, c_min, delta, is_somatic);
    return ret.posterior_error;
}

//
// Main
//
int graphConcordanceMain(int argc, char** argv)
{
    parseGraphConcordanceOptions(argc, argv);

    Timer* pTimer = new Timer(PROGRAM_IDENT);

    // Load the reference
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Index the germline variants so we can make germline haplotypes
    std::cerr << "Loading germline variant index..." << std::flush;
    VariantIndex germlineIndex(opt::germlineFile, refTable);
    std::cerr << "done\n";

    // Load FM-index of the reads
    BWTIndexSet variantIndex;
    std::cerr << "Loading variant read index... " << std::flush;
    std::string variantPrefix = stripGzippedExtension(opt::variantFile);
    variantIndex.pBWT = new BWT(variantPrefix + BWT_EXT, 256);
    variantIndex.pSSA = new SampledSuffixArray(variantPrefix + SAI_EXT, SSA_FT_SAI);
    variantIndex.pCache = new BWTIntervalCache(11, variantIndex.pBWT);
    variantIndex.pQualityTable = new QualityTable;
    std::cerr << "done\n";
    
    //
    BWTIndexSet baseIndex;
    std::cerr << "Loading base read index index... " << std::flush;
    std::string basePrefix = stripGzippedExtension(opt::baseFile);
    baseIndex.pBWT = new BWT(basePrefix + BWT_EXT, 256);
    baseIndex.pCache = new BWTIntervalCache(11, baseIndex.pBWT);
    std::cerr << "done\n";

    // Calculate the parameters for the k-mer model
    ModelParameters variant_params = calculate_model_parameters(opt::k, 50000, variantIndex);
    GenomeEstimates variant_genome = estimate_genome_size_from_k_counts(opt::k, variantIndex);
    
    fprintf(stderr, "Variant -- Kmer mode: %lf\n", variant_params.mode);
    fprintf(stderr, "Variant -- Genome size: %zu\n", variant_genome.genome_size);
    
    ModelParameters base_params = calculate_model_parameters(opt::k, 50000, baseIndex);
    GenomeEstimates base_genome = estimate_genome_size_from_k_counts(opt::k, baseIndex);
    
    fprintf(stderr, "Base -- Kmer mode: %lf\n", base_params.mode);
    fprintf(stderr, "Base -- Genome size: %zu\n", base_genome.genome_size);

    std::ifstream input(opt::vcfFile.c_str());
    std::string line;
    while(getline(input, line))
    {
        if(line.empty())
            continue;

        if(line[0] == '#')
        {
            std::cout << line << "\n";
            continue;
        }
        
        // parse record
        VCFRecord record(line);
        if(record.isMultiAllelic())
        {
            std::cerr << "Error: multi-allelic VCF found, please run vcfbreakmulti\n";
            exit(EXIT_FAILURE);
        }

        // Grab nearby variants
        size_t flanking_size = opt::k;
        VariantRecordVector nearby_vector = germlineIndex.getNearVariants(record.refName, 
                                                                          record.refPosition, 
                                                                          flanking_size);
        
        if(opt::verbose > 0)
        {
            std::cerr << "\n\n Variant: " << record << "\n";
            std::cerr << "Nearby: " << nearby_vector.size() << "\n";
        }

        std::string classification = "UNKNOWN";

        size_t local_k = flanking_size;
        ReconstructedHaplotype rc;
        if( nearby_vector.size() > 4 || 
            pickLocalKAndReconstruct(record, nearby_vector, refTable, variantIndex, baseIndex, flanking_size,
                                    variant_params, base_params, rc, local_k) != 0 )
        {
            record.setQuality(0);
            classification = "TOO_COMPLEX";
        }
        else
        {
            std::string left_x = rc.left_base_kmer;
            std::string left_y = rc.left_somatic_kmer;
            std::string left_z = get_neighbor(left_y, ED_SENSE, variantIndex);

            std::string right_x = rc.right_base_kmer;
            std::string right_y = rc.right_somatic_kmer;
            std::string right_z = get_neighbor(right_y, ED_ANTISENSE, variantIndex);

            int sum_y = kmerCount(left_y, variantIndex) + kmerCount(right_y, variantIndex);
            int sum_z = kmerCount(left_z, variantIndex) + kmerCount(right_z, variantIndex);

            // Generate various output statistics
            
            record.addComment("OrigQual", record.quality);
            record.addComment("HapUniqueKmers", (int)rc.unique_kmers);
            record.addComment("LocalK", (int)local_k);

            //
            // Number of ALT reads from freebayes
            //
            if(!record.sampleStr.empty())
            {
                StringVector gt_fields = split(record.sampleStr[0], ':');
                if(gt_fields.size() >= 6)
                {
                    std::stringstream parser(gt_fields[5]);
                    int freebayes_alt_reads = 0;
                    parser >> freebayes_alt_reads;
                    record.addComment("MappedAlt", freebayes_alt_reads);
                    
                    int freebayes_total_reads = 0;
                    std::stringstream parser2(gt_fields[2]);
                    parser2 >> freebayes_total_reads;
                    record.addComment("MappedTotal", freebayes_total_reads);
                }
            }

            //
            // Probability model
            //
            double p_left_error_variant = 0.99;
            double p_right_error_variant = 0.99;
            
            double p_left_error_base = 0.99;
            double p_right_error_base = 0.99;

            if(!left_x.empty() && !left_y.empty())
            {
                p_left_error_variant = apply_branch_model(left_x, left_y, left_z,
                                                          variantIndex, variant_params, variant_genome, true);

                p_left_error_base = apply_branch_model(left_x, left_y, left_z,
                                                       baseIndex, base_params, base_genome, false);
            }

            if(!right_x.empty() && !right_y.empty())
            {
                p_right_error_variant = apply_branch_model(right_x, right_y, right_z,
                                                           variantIndex, variant_params, variant_genome, true);

                p_right_error_base = apply_branch_model(right_x, right_y, right_z,
                                                        baseIndex, base_params, base_genome, false);
            }
            
            double p_error_variant = p_left_error_variant * p_right_error_variant;
            double p_error_base = p_left_error_base * p_right_error_base;
            double p_somatic = (1.0 - p_error_variant) * p_error_base;

            double model_qual = -10.0 * log(p_error_variant) / log(10);
            record.addComment("ModelQual", model_qual);

            double prior_error = 0.01;
            double log_p_error_variant = SGAStats::logBinomial(sum_y, sum_y + sum_z, prior_error);
            double binomial_qual = -10.0 * log_p_error_variant / log(10);
            record.addComment("KmerAlt", sum_y);
            record.addComment("KmerTotal", sum_y + sum_z);
            record.addComment("BinomialQual", binomial_qual);

            // QC
            int sum_y_ss = kmerCount(left_y, variantIndex) + 
                           kmerCount(right_y, variantIndex);
            bool is_balanced = sum_y_ss < sum_y && sum_y_ss > 0;

            //
            // Set VCF Qual
            //
            double qual = 0.0;
            
            int left_y_base = kmerCount(left_y, baseIndex);
            int right_y_base = kmerCount(right_y, baseIndex);
            
            int left_y_variant  = kmerCount(left_y, variantIndex);
            int right_y_variant = kmerCount(right_y, variantIndex);

            record.addComment("YLeftBase", left_y_base);
            record.addComment("YRightBase", right_y_base);
            
            record.addComment("YLeftVariant", left_y_variant);
            record.addComment("YRightVariant", right_y_variant);

            // Homopolymer length
            int hplen = calculateHomopolymerLength(record, &refTable);
            record.addComment("HPLen", hplen);

            // Dust
            double dust = HapgenUtil::calculateDustScoreAtPosition(record.refName, 
                                                                   record.refPosition, 
                                                                   &refTable);

            int sum_y_base = left_y_base + right_y_base;
            // Note: here we're allowing _1_ sum_y_base.  This shouldn't be a
            // hardcoded number.
            if(sum_y_base < 2 && hplen <= 7 && dust <= 2.0)
                qual = model_qual;

            //double qual_base_not_error = -10.0 * log(1 - p_error_base) / log(10);
            //double qual = std::min(qual_variant, qual_base_not_error);

            fprintf(stderr, "p_e_v: %.3lf\n", exp(log_p_error_variant));
            //fprintf(stderr, "p_e_b: %.3lf\n", p_error_base);
            //fprintf(stderr, "p_somatic: %.3lf\n", p_somatic);


            //double qual = -10.0 * log(1 - p_somatic) / log(10);
            record.setQuality(qual);
            fprintf(stderr, "Qual: %.1lf\n", qual);
        }
        record.addComment("KmerClassification", classification);
        std::cout << record << "\n";
    }
    
    // Cleanup
    delete variantIndex.pBWT;
    delete variantIndex.pCache;
    delete baseIndex.pBWT;
    delete baseIndex.pCache;
    delete pTimer;

    return 0;
}

// 
// Handle command line arguments
//
void parseGraphConcordanceOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'b': arg >> opt::baseFile; break;
            case 'r': arg >> opt::variantFile; break;
            case 'g': arg >> opt::germlineFile; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GRAPH_CONCORDANCE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GRAPH_CONCORDANCE_VERSION_MESSAGE;
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

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": a reference file must be given.\n";
        die = true;
    }

    if(opt::germlineFile.empty())
    {
        std::cerr << SUBPROGRAM ": a germline variants file must be given.\n";
        die = true;
    }

    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << GRAPH_CONCORDANCE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
