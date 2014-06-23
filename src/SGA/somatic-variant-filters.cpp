//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// variant-filtering - apply various filters
// to a somatic VCF file
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "somatic-variant-filters.h"
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
#include "api/BamReader.h"

// Types
typedef HashMap<std::string, std::string> StringStringHash;

//
// Getopt
//
#define SUBPROGRAM "somatic-variant-filters"
static const char *SOMATIC_VARIANT_FILTERS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2014 Ontario Institute for Cancer Research\n";

static const char *SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... VCF_FILE\n"
"Count read support for variants in a vcf file\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"          --reference=STR              load the reference genome from FILE\n"
"          --tumor-bam=STR              load the aligned tumor reads from FILE\n"
"          --normal-bam=STR             load the aligned normal reads from FILE\n"
"\n"
"Filtering cutoffs:\n"
"          --min-af=FLOAT               minimum allele frequency (AF tag)\n"
"          --min-var-dp=INT             minimum variant depth (VarDP)\n"
"          --max-strand-bias=FLOAT      maximum strand bias (SB)\n"
"          --max-hp-length=INT          maximum homopolymer context length (HPlen)\n"
"          --max-dust=FLOAT             maximum dust scores (Dust)\n"
"          --max-normal-reads=INT       maximum normal reads showing variant\n"
"          --min-normal-depth=INT       minimum normal depth at the locus\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string vcfFile;
    static std::string referenceFile;
    static std::string tumorBamFile;
    static std::string normalBamFile;

    static double minAF = 0.0f;
    static int minVarDP = 0;
    static double maxStrandBias = 2.0f;
    static int maxHPLen = 7;
    static double maxDust = 2.0f;
    static size_t maxNormalReads = 1;
    static size_t minNormalDepth = 5;
    static int minMedianQuality = 15;
}

static const char* shortopts = "t:v";

enum { OPT_HELP = 1, 
       OPT_VERSION, 
       OPT_REFERENCE, 
       OPT_TUMOR_BAM, 
       OPT_NORMAL_BAM,
       OPT_MIN_AF,
       OPT_MIN_VARDP,
       OPT_MAX_SB,
       OPT_MAX_HP,
       OPT_MAX_DUST,
       OPT_MAX_NORMAL_READS,
       OPT_MIN_NORMAL_DEPTH };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "reference",             required_argument, NULL, OPT_REFERENCE },
    { "tumor-bam",             required_argument, NULL, OPT_TUMOR_BAM },
    { "normal-bam",            required_argument, NULL, OPT_NORMAL_BAM },
    { "min-af",                required_argument, NULL, OPT_MIN_AF },
    { "min-var-dp",            required_argument, NULL, OPT_MIN_VARDP },
    { "max-strand-bias",       required_argument, NULL, OPT_MAX_SB },
    { "max-hp-length",         required_argument, NULL, OPT_MAX_HP },
    { "max-dust",              required_argument, NULL, OPT_MAX_DUST },
    { "max-normal-reads",      required_argument, NULL, OPT_MAX_NORMAL_READS },
    { "min-normal-depth",      required_argument, NULL, OPT_MIN_NORMAL_DEPTH },
    { "help",                  required_argument, NULL, OPT_HELP },
    { "version",               required_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

struct CoverageStats
{
    CoverageStats() : n_total_reads(0), n_evidence_reads(0) {} 
    size_t n_total_reads;
    size_t n_evidence_reads;
    std::vector<int> snv_evidence_quals;
};

// This struct holds segments
// of a read that match a variant
// The read is broken into (up to) 3 pieces
//  -the segment that aligns to the reference before the variant
//  -the segment that aligns to the variant region
//  -the segment that aligns after the variant
struct VariantReadSegments
{
    std::string preSegment;
    std::string variantSegment;
    std::string postSegment;

    std::string preQual;
    std::string variantQual;
    std::string postQual;
};

VariantReadSegments splitReadAtVariant(const BamTools::BamAlignment& alignment,
                                       const VCFRecord& record)
{
    // Set up a vector of aligned reference positions
    // for every base of the read
    std::vector<int> reference_positions(alignment.QueryBases.size(), -1);
    size_t read_idx = 0;
    size_t ref_idx = alignment.Position;

    std::stringstream cigar;
    // 
    for(size_t ci = 0; ci < alignment.CigarData.size(); ++ci)
    {
        BamTools::CigarOp op = alignment.CigarData[ci];
        switch(op.Type) {
            case 'M':
                for(size_t l = 0; l < op.Length; ++l)
                    reference_positions[read_idx++] = ref_idx++;
                break;
            case 'D':
            case 'N':
                ref_idx += op.Length;
                break;
            case 'I':
            case 'S':
                read_idx += op.Length;
                break;
        }

        cigar << op.Length << op.Type;
    }

    // Find the last aligned position before the variant and the first position after
    int rightmost_before_variant = -1;
    int leftmost_after_variant = std::numeric_limits<int>::max();

    int variant_start = record.refPosition - 1; // VCF is 1-based, bamtools is 0-based
    int variant_end = variant_start + record.refStr.size();

    for(int i = 0; i < (int)reference_positions.size(); ++i) {
        if(reference_positions[i] != -1 && reference_positions[i] < variant_start && i > rightmost_before_variant)
            rightmost_before_variant = i;
        if(reference_positions[i] != -1 && reference_positions[i] >= variant_end && i < leftmost_after_variant)
            leftmost_after_variant = i;
    }
    //fprintf(stderr, "rightmost: %d leftmost_after_variant: %d\n", rightmost_before_variant, leftmost_after_variant);

    VariantReadSegments ret;
    ret.preSegment = alignment.QueryBases.substr(0, rightmost_before_variant + 1);
    ret.preQual = alignment.Qualities.substr(0, rightmost_before_variant + 1);

    ret.variantSegment = alignment.QueryBases.substr(rightmost_before_variant + 1, leftmost_after_variant - rightmost_before_variant - 1);
    ret.variantQual = alignment.Qualities.substr(rightmost_before_variant + 1, leftmost_after_variant - rightmost_before_variant - 1);

    if(leftmost_after_variant < (int)alignment.QueryBases.size())
    {
        ret.postSegment = alignment.QueryBases.substr(leftmost_after_variant);
        ret.postQual = alignment.Qualities.substr(leftmost_after_variant);
    }

    return ret;
}

//
//
//
CoverageStats getVariantCoverage(BamTools::BamReader* pReader, const VCFRecord& record, const ReadTable* refTable)
{
    CoverageStats stats;
    
    static const int flankingSize = 100;
    static const double minPercentIdentity = 95.0f;

    bool is_snv = record.refStr.size() == 1 && record.varStr.size() == 1;

    // Grab the reference haplotype
    int eventLength = record.varStr.length();
    int zeroBasedPos = record.refPosition - 1;
    int start = zeroBasedPos - flankingSize - 1;
    if(start < 0)
        start = 0;

    int end = zeroBasedPos + eventLength + 2 * flankingSize;
    const SeqItem& chr = refTable->getRead(record.refName);
    if(end > (int)chr.seq.length())
        end = (int)chr.seq.length();

    std::string reference_haplotype = chr.seq.substr(start, end - start);
    int translatedPos = zeroBasedPos - start;

    std::string variant_haplotype = reference_haplotype;
    
    // Ensure that the reference string at the variant matches the expected
    assert(variant_haplotype.substr(translatedPos, record.refStr.length()) == record.refStr);
    variant_haplotype.replace(translatedPos, record.refStr.length(), record.varStr);

    // Grab all reads in reference region
    int refID = pReader->GetReferenceID(record.refName);
    if(refID < 0)
        return stats;

    int refStart = record.refPosition;
    int refEnd = record.refPosition;
    pReader->SetRegion(refID, refStart, refID, refEnd);
    BamTools::BamAlignment alignment;

    while(pReader->GetNextAlignment(alignment)) {
        
        VariantReadSegments segments = splitReadAtVariant(alignment, record);

        if(opt::verbose > 1)
        {
            fprintf(stderr, "var: %zu %s -> %s\n",  record.refPosition, record.refStr.c_str(), record.varStr.c_str());
            fprintf(stderr, "pos: %d\n",  alignment.Position);
            fprintf(stderr, "strand: %s\n", alignment.IsReverseStrand() ? "-" : "+");
            fprintf(stderr, "read: %s\n", alignment.QueryBases.c_str());
            fprintf(stderr, "qual: %s\n", alignment.Qualities.c_str());
            fprintf(stderr, "alnb: %s\n", alignment.AlignedBases.c_str());
            
            fprintf(stderr, "Pre: %s\n",  segments.preSegment.c_str());
            fprintf(stderr, "Var: %s\n",  segments.variantSegment.c_str());
            fprintf(stderr, "Pos: %s\n",  segments.postSegment.c_str());
            
            fprintf(stderr, "PreQual: %s\n",  segments.preQual.c_str());
            fprintf(stderr, "VarQual: %s\n",  segments.variantQual.c_str());
            fprintf(stderr, "PosQual: %s\n",  segments.postQual.c_str());
        }

        if( (segments.preSegment.size() > 0 || segments.postSegment.size() > 0) &&
            segments.variantSegment.size() > 0)
        {
            stats.n_total_reads += 1;
            
            // Align the read to the reference and variant haplotype
            SequenceOverlap ref_overlap = Overlapper::computeOverlapAffine(alignment.QueryBases, reference_haplotype);
            SequenceOverlap var_overlap = Overlapper::computeOverlapAffine(alignment.QueryBases, variant_haplotype);
            
            /*
            std::cout << "OverlapRef: \n";
            ref_overlap.printAlignment(alignment.QueryBases, reference_haplotype);
            std::cout << "OverlapVar: \n"; 
            var_overlap.printAlignment(alignment.QueryBases, variant_haplotype);
            */
            
            bool quality_alignment = (ref_overlap.getPercentIdentity() >= minPercentIdentity || 
                                     var_overlap.getPercentIdentity() >= minPercentIdentity);

            bool is_evidence_read = quality_alignment && var_overlap.score > ref_overlap.score;
            if(is_evidence_read)
            {
                stats.n_evidence_reads += 1;
                if(is_snv)
                {
                    if(segments.variantQual.size() == 1)
                    {
                        char qb = segments.variantQual[0];
                        int q = Quality::char2phred(qb);
                        stats.snv_evidence_quals.push_back(q);
                    }
                }
            }
        }
    }

    return stats;
}

void makeTagHash(const VCFRecord& record, StringStringHash& tagHash)
{
    for(size_t i = 0; i < record.comments.size(); i++)
    {
        StringVector kv = split(record.comments[i], '=');
        assert(kv.size() == 2);
        tagHash[kv[0]] = kv[1];
    }
}

template<typename T>
bool getTagValue(StringStringHash& tagHash, const std::string& key, T& out)
{
    StringStringHash::iterator iter = tagHash.find(key);
    if(iter == tagHash.end())
        return false; // tag not available
    std::stringstream parser(iter->second);
    parser >> out;
    return true;
}

std::string makeVariantKey(const VCFRecord& record)
{
    std::stringstream ss;
    ss << record.refName << ":" << record.refPosition << ":"
       << record.refStr << ":" << record.varStr;
    return ss.str(); 
}

//
// Main
//
int somaticVariantFiltersMain(int argc, char** argv)
{
    parseSomaticVariantFiltersOptions(argc, argv);

    Timer* pTimer = new Timer(PROGRAM_IDENT);
    
    // Load Reference
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Load BAMs
    BamTools::BamReader* pTumorBamReader = new BamTools::BamReader;
    pTumorBamReader->Open(opt::tumorBamFile);
    pTumorBamReader->LocateIndex();

    assert(pTumorBamReader->HasIndex());

    BamTools::BamReader* pNormalBamReader = new BamTools::BamReader;
    pNormalBamReader->Open(opt::normalBamFile);
    pNormalBamReader->LocateIndex();
    assert(pNormalBamReader->HasIndex());

    // Track duplicated variants
    HashSet<std::string> duplicateHash;

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

        // Check if we've seen this variant already
        std::string key = makeVariantKey(record);
        if(duplicateHash.find(key) != duplicateHash.end())
            continue;
        else
            duplicateHash.insert(key);

        if(opt::verbose > 0)
        {
            std::stringstream ss;
            ss << "Variant: " << record << "\n";
            fprintf(stderr, "===============================================\n%s", ss.str().c_str());
        }

        StringStringHash tagHash;
        makeTagHash(record, tagHash);

        StringVector fail_reasons;

        int hplen;
        if(getTagValue(tagHash, "HPLen", hplen) && hplen > opt::maxHPLen)
            fail_reasons.push_back("Homopolymer");

        double dust;
        if(getTagValue(tagHash, "Dust", dust) && dust > opt::maxDust)
            fail_reasons.push_back("LowComplexity");
        
        double af;
        if(getTagValue(tagHash, "AF", af) && af < opt::minAF)
            fail_reasons.push_back("LowAlleleFrequency");

        int varDP;
        if(getTagValue(tagHash, "VarDP", varDP) && varDP < opt::minVarDP)
            fail_reasons.push_back("LowVarDP");

        double strandBias;
        if(getTagValue(tagHash, "SB", strandBias) && strandBias >= opt::maxStrandBias)
            fail_reasons.push_back("StrandBias");

        CoverageStats tumor_stats = getVariantCoverage(pTumorBamReader, record, &refTable);
        CoverageStats normal_stats = getVariantCoverage(pNormalBamReader, record, &refTable);

        if(opt::verbose > 0)
        {
            fprintf(stderr, "Tumor: [%zu %zu]\n",  tumor_stats.n_total_reads, tumor_stats.n_evidence_reads);
            fprintf(stderr, "Normal: [%zu %zu]\n", normal_stats.n_total_reads, normal_stats.n_evidence_reads);
        }

        if(normal_stats.n_evidence_reads > opt::maxNormalReads)
            fail_reasons.push_back("NormalEvidence");
        
        if(normal_stats.n_total_reads < opt::minNormalDepth)
            fail_reasons.push_back("LowNormalDepth");

        if(!tumor_stats.snv_evidence_quals.empty())
        {
            std::sort(tumor_stats.snv_evidence_quals.begin(), tumor_stats.snv_evidence_quals.end());

            double median_quality = 0.0f;
            if(tumor_stats.snv_evidence_quals.size() % 2 == 1)
            {
                int m = tumor_stats.snv_evidence_quals.size() / 2;
                median_quality = tumor_stats.snv_evidence_quals[m];
            }
            else
            {
                int m = tumor_stats.snv_evidence_quals.size() / 2;
                median_quality = (tumor_stats.snv_evidence_quals[m - 1] + tumor_stats.snv_evidence_quals[m]) / 2.0f;
            }
        
            if(median_quality < opt::minMedianQuality)
                fail_reasons.push_back("LowQuality");
        }

        if(!fail_reasons.empty())
        {
            if(record.passStr != "PASS" || record.passStr != ".")
                fail_reasons.insert(fail_reasons.begin(), record.passStr);

            std::stringstream strss;
            std::copy(fail_reasons.begin(), fail_reasons.end(), std::ostream_iterator<std::string>(strss, ";"));
            record.passStr = strss.str();
            record.passStr.erase(record.passStr.size() - 1); // erase trailing ;
        }

        std::cout << record << "\n";
    }
    
    // Cleanup
    delete pTumorBamReader;
    delete pNormalBamReader;
    delete pTimer;

    return 0;
}

// 
// Handle command line arguments
//
void parseSomaticVariantFiltersOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_TUMOR_BAM: arg >> opt::tumorBamFile; break;
            case OPT_NORMAL_BAM: arg >> opt::normalBamFile; break;
            case OPT_MIN_AF: arg >> opt::minAF; break;
            case OPT_MIN_VARDP: arg >> opt::minVarDP; break;
            case OPT_MAX_SB: arg >> opt::maxStrandBias; break;
            case OPT_MAX_HP: arg >> opt::maxHPLen; break;
            case OPT_MAX_DUST: arg >> opt::maxDust; break;
            case OPT_MAX_NORMAL_READS: arg >> opt::maxNormalReads; break;
            case OPT_MIN_NORMAL_DEPTH: arg >> opt::minNormalDepth; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SOMATIC_VARIANT_FILTERS_VERSION_MESSAGE;
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

    if(opt::tumorBamFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --tumor-bam must be provided\n";
        die = true;
    }
    
    if(opt::normalBamFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --normal-bam must be provided\n";
        die = true;
    }

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --reference must be provided\n";
        die = true;
    }

    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
