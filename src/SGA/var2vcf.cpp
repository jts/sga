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
void parseVariants(const ReadTable* pRefTable, const BamTools::BamReader* pReader, const BamRecordVector& records);

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
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static std::string referenceFile;
    static std::string bamFile;
}

static const char* shortopts = "r:ov";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "refrence",        required_argument, NULL, 'r' },
    { "outfile",         required_argument, NULL, 'o' },
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

    Timer* pTimer = new Timer(PROGRAM_IDENT);    
    
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

    std::cout << "Reading from BAM: " << opt::bamFile << "\n";
    std::cout << "Reading from ref: " << opt::referenceFile << "\n";

    VariantGroupReader groupReader(pBamReader);
    while(1)
    {
        BamRecordVector records;
        if(groupReader.readVariantGroup(records))
        {
            parseVariants(&refTable, pBamReader, records);
        }
        else
        {
            break;
        }
    }

    pBamReader->Close();
    delete pTimer;
    delete pBamReader;
    return 0;
}

void parseVariants(const ReadTable* pRefTable, const BamTools::BamReader* pReader, const BamRecordVector& records)
{
    // classify the sequences in the vector as base 
    // or variant
    std::string baseString;
    std::string varString;

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

    printf("Ref: %s coords: [%zu %zu] Name: %s\n", refName.c_str(), refStartPos, refEndPos, records[0].Name.c_str());
    printf("base: %s\n", baseString.c_str());
    printf(" var: %s\n", varString.c_str());
    printf(" ref: %s\n", refString.c_str());

    std::cout << "Alignments\n";
    StdAlnTools::globalAlignment(refString, baseString, true);
    StdAlnTools::globalAlignment(refString, varString, true);
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
