//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// var2vcf - convert variants from sga graph-diff
// or sga assembly to a vcf file of differences
// with respect to a reference
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include "ReadTable.h"
#include "Util.h"
#include "var2vcf.h"
#include "Timer.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "StdAlnTools.h"
#include "VCFUtil.h"

//
typedef std::vector<BamTools::BamAlignment> BamRecordVector;
typedef std::map<std::string, std::string> StrStrMap;

// Small class to read alignments of a variant group from a bam file
struct VariantGroupReader
{
    public:
        
        VariantGroupReader(BamTools::BamReader* pReader) : m_pReader(pReader)
        {
            // prime the reader with a single record
            bool success = m_pReader->GetNextAlignment(m_cachedAlignment);
            if(!success)
            {
                std::cerr << "Error: no records in the BAM file\n";
                exit(EXIT_FAILURE);
            }
            m_doneRead = false;
        }

        // Read a group of records from the bam belonging
        // to the same variant. Returns true if records
        // is successfully populated with alignments
        bool readVariantGroup(BamRecordVector& records)
        {
            if(m_doneRead)
                return false;

            std::string cachedName = m_cachedAlignment.Name;
            size_t sufPos = cachedName.find_last_of('-');
            assert(sufPos != std::string::npos && sufPos > 0);
            std::string targetSuffix = cachedName.substr(sufPos);

            // Read alignments from the string as long as they have the same
            // suffix as the cached alignment
            records.push_back(m_cachedAlignment);
            while(1)
            {
                if(!m_pReader->GetNextAlignment(m_cachedAlignment))
                {
                    m_doneRead = true;
                    break; // no more reads but the last group should be processed
                }

                if(m_cachedAlignment.Name.rfind(targetSuffix) != std::string::npos)
                    records.push_back(m_cachedAlignment); // part of the group
                else
                    break; // the read parsed is not part of the same group
            } 
            return true;
        }

    private:
        bool m_doneRead;
        BamTools::BamReader* m_pReader;
        BamTools::BamAlignment m_cachedAlignment;
};


// Functions
VCFReturnCode preprocessVariants(BamRecordVector& records);

VCFReturnCode parseVariants(const ReadTable* pRefTable, const BamTools::BamReader* pReader, 
                            const BamRecordVector& records, VCFVector& outVCFRecords);

//
// Getopt
//
#define SUBPROGRAM "var2vcf"
static const char *VAR2VCF_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *VAR2VCF_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... BAMFILE\n"
"Convert the sequence aligned to a reference in BAMFILE into a VCF file of differences\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reference=FILE             read the reference sequences from FILE\n"
"      -o, --outfile=FILE               write the results to FILE\n"
"      -q, --min-quality=Q              discard variants with mapping quality less than Q (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static std::string referenceFile;
    static std::string bamFile;
    static int minQuality = 1;
    static int exactMatchRequired = 21;
}

static const char* shortopts = "r:o:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "refrence",        required_argument, NULL, 'r' },
    { "outfile",         required_argument, NULL, 'o' },
    { "min-quality",     required_argument, NULL, 'q' },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int var2vcfMain(int argc, char** argv)
{
    parseVar2VCFOptions(argc, argv);

    // Read the reference
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Open the bam files for reading/writing
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    if(!pBamReader->Open(opt::bamFile))
    {
        std::cerr << "Failed to open BAM file: " << opt::bamFile << "\n";
        exit(EXIT_FAILURE);
    }
    WARN_ONCE("Add dust check");
    //
    // Read the alignments from the BAM file as a group of variants
    // and convert them to VCF
    //
    VariantGroupReader groupReader(pBamReader);
    IntVector returnCodeStats(VCF_NUM_RETURN_CODES, 0);
 
    VCFVector vcfRecords;
    size_t numGroups = 0;
    while(1)
    {
        BamRecordVector records;

        bool groupOK = groupReader.readVariantGroup(records);
        if(!groupOK)
            break;

        numGroups += 1;
        VCFReturnCode code = preprocessVariants(records);
        if(code != VCF_OK)
        {
            // Failed preprocess checks
            returnCodeStats[code] += 1;
            continue;
        }
        
        code = parseVariants(&refTable, pBamReader, records, vcfRecords);
        returnCodeStats[code] += 1;
    }

    // Sort the VCF records by chromosome then position
    std::sort(vcfRecords.begin(), vcfRecords.end(), VCFRecord::sort);
    
    //
    // Output VCF
    //

    // Build program string
    std::stringstream progSS;
    progSS << "sga";
    for(int i = 0; i < argc; ++i)
        progSS << " " << argv[i];

    IntVector classificationStats(VCF_NUM_CLASSIFICATIONS, 0);
    std::ostream* pWriter = createWriter(opt::outFile);
    VCFUtil::writeHeader(pWriter, PROGRAM_IDENT, stripDirectories(opt::bamFile), stripDirectories(opt::referenceFile));
    for(size_t i = 0; i < vcfRecords.size(); ++i)
    {
        VCFClassification code = vcfRecords[i].classify();
        assert(code < (int)classificationStats.size());
        classificationStats[code] += 1;
        *pWriter << vcfRecords[i] << "\n";
    }
    delete pWriter;

    // Print stats
    printf("Total variant sequences: %zu\n", numGroups);
    printf(" -- Successfully converted: %d\n", returnCodeStats[VCF_OK]);
    printf(" -- Failed due to flanking exact match check: %d\n", returnCodeStats[VCF_EXACT_MATCH_FAILED]);
    printf(" -- Failed due to low mapping quality: %d\n", returnCodeStats[VCF_MAP_QUALITY_FAILED]);
    printf(" -- Failed due to ambiguous base sequence mappping: %d\n", returnCodeStats[VCF_BASE_MULTIMAP_FAILED]);
    printf(" -- Failed due to partially mapped base sequence: %d\n", returnCodeStats[VCF_BASE_PARTIALMAP_FAILED]);
    printf(" -- Failed due to invalid multiple alignment: %d\n", returnCodeStats[VCF_INVALID_MULTIALIGNMENT]);
    printf("\n");
    printf("Totat variants: %zu\n", vcfRecords.size());
    printf(" -- substitutions: %d\n", classificationStats[VCF_SUB]);
    printf(" -- deletions: %d\n", classificationStats[VCF_DEL]);
    printf(" -- insertions: %d\n", classificationStats[VCF_INS]);
    printf(" -- complex: %d\n", classificationStats[VCF_COMPLEX]);
    pBamReader->Close();
    delete pBamReader;
    return 0;
}

// Perform sanity checks and filtering of the records
VCFReturnCode preprocessVariants(BamRecordVector& records)
{
    // Perform quality check   
    for(size_t i = 0; i < records.size(); ++i)
    {
        // Only check the base sequence's map quality
        // The variant can potentially be so different from the reference
        // that it does not map
        if(records[i].Name.find("base") != std::string::npos)
        {
            if(records[i].MapQuality < opt::minQuality)
                return VCF_MAP_QUALITY_FAILED;

            // Check if the alignment is softclipped
            assert(!records[i].CigarData.empty());
            if(records[i].CigarData.front().Type == 'S' || 
               records[i].CigarData.back().Type == 'S')
            {
                return VCF_BASE_PARTIALMAP_FAILED;
            }
        }
    }

    // Perform duplicate mapping check
    int numBaseRecords = 0;
    int numVariantRecords = 0;
    for(size_t i = 0; i < records.size(); ++i)
    {
        if(records[i].Name.find("base") == 0)
            numBaseRecords += 1;
        else if(records[i].Name.find("variant") == 0)
            numVariantRecords += 1;
        else
        {
            std::cerr << "Error: Unexpected record name: " << records[i].Name << "\n";
            exit(EXIT_FAILURE);
        }
    }
    
    // Fail if the base is multimapped
    if(numBaseRecords > 1)
        return VCF_BASE_MULTIMAP_FAILED;

    // Map a new copy of the records with only 1 entry for the variant
    if(numVariantRecords > 1)
    {
        BamRecordVector outRecords;
        // Make a copy of the records with only 1 variant record
        for(size_t i = 0; i < records.size(); ++i)
        {
            if(records[i].Name.find("base") == 0)
            {
                outRecords.push_back(records[i]);
            }
            else if(records[i].Name.find("variant") == 0)
            {
                outRecords.push_back(records[i]);
                break;
            }
        }
        records = outRecords;
    }

    return VCF_OK;
}

// Perform the actual conversion of the collection of BAM records into a call
VCFReturnCode parseVariants(const ReadTable* pRefTable, const BamTools::BamReader* pReader, 
                            const BamRecordVector& records, VCFVector& outVCFRecords)
{
    // classify the sequences in the vector as base 
    // or variant
    std::string baseString;
    std::string varString;
    std::string varName;
    std::string refName;
    size_t refStartPos = 0;
    size_t refEndPos = 0;
    bool bBaseIsRC = false;

    if(records.size() != 2)
    {
        std::cerr << "Error, expected 2 records (got: " << records.size() << ")\n";
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < records.size(); ++i)
    {
        if(records[i].Name.find("base") != std::string::npos)
        {
            if(!refName.empty())
            {
                std::cerr << "Error: base sequence " << records[i].Name << " does not have a single alignment\n";
                exit(EXIT_FAILURE);
            }

            assert(records[i].IsMapped());
            baseString = records[i].QueryBases;
            refStartPos = records[i].Position;
            refEndPos = records[i].GetEndPosition();
            const BamTools::RefVector& refVector = pReader->GetReferenceData();
            assert(records[i].RefID < (int)refVector.size());
            refName = refVector[records[i].RefID].RefName;
            bBaseIsRC = records[i].IsReverseStrand();
        }

        if(records[i].Name.find("variant") != std::string::npos)
        {
            varString = records[i].QueryBases;
            varName = records[i].Name;

            // We always keep the variant string in its original
            // strand, then flip it to match the base
            if(records[i].IsReverseStrand())
                varString = reverseComplement(varString);
        }
    }

    if(varString.empty() || baseString.empty())
    {
        std::cerr << "Error variant string or base string not found\n";
        exit(EXIT_FAILURE);
    }

    // Correct the orientation of the var string to match the base
    if(bBaseIsRC)
        varString = reverseComplement(varString);

    // Extract the portion of the reference
    const SeqItem& refItem = pRefTable->getRead(refName);
    assert(refStartPos < refItem.seq.length() && refEndPos < refItem.seq.length());
    std::string refString = refItem.seq.substr(refStartPos, refEndPos - refStartPos + 1);
    
    // Perform the actual conversion
    VCFReturnCode code = VCFUtil::generateVCFFromCancerVariant(refString, baseString, varString, refName, refStartPos, varName, opt::exactMatchRequired, opt::verbose, outVCFRecords);
    return code;
}

// 
// Handle command line arguments
//
void parseVar2VCFOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'r': arg >> opt::referenceFile; break;
            case 'q': arg >> opt::minQuality; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << VAR2VCF_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << VAR2VCF_VERSION_MESSAGE;
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
        std::cerr << SUBPROGRAM ": a reference file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << VAR2VCF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::bamFile = argv[optind++];

    if(opt::outFile.empty())
        opt::outFile = "variants.vcf";
}
