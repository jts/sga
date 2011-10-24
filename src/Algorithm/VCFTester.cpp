///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFTester - Run dindel on each variant in a vcf file
//
#include "VCFTester.h"
#include "StdAlnTools.h"
#include "DindelUtil.h"
//
//
//
VCFTester::VCFTester(const GraphCompareParameters& params) : m_parameters(params), 
                                                             m_baseVCFFile("vtest.base.vcf","w"),
                                                             m_variantVCFFile("vtest.variant.vcf","w")
{
    m_baseVCFFile.outputHeader("stub", "stub");
    m_variantVCFFile.outputHeader("stub", "stub");
    DindelUtil::initializeCodeCounts(m_returnCodes);
}

//
VCFTester::~VCFTester()
{
    std::cout << "Done testing variants\n";
    DindelUtil::printReturnReport(m_returnCodes);
}

void VCFTester::process(const VCFFile::VCFEntry& record)
{
    // Extract the reference string for this record
    // We want a full kmer to the left and right of the variant
    int eventLength = record.alt.length();

    int zeroBasedPos = record.pos - 1;
    int start = zeroBasedPos - m_parameters.kmer - 1;
    int end = zeroBasedPos + eventLength + m_parameters.kmer;

    const SeqItem& chr = m_parameters.pRefTable->getRead(record.chrom);

    std::string refStr = chr.seq.substr(start, end - start);

    // Apply the variant to the reference string
    // Update the coordinate to be wrt the start of the reference substring
    int relative_pos = zeroBasedPos - start;
    std::string varStr = applyVariant(refStr, relative_pos, record.ref, record.alt);

    std::cout << "Testing variant";
    record.write(std::cout);

    StdAlnTools::globalAlignment(refStr, varStr, true);

    // Run dindel
    DindelReturnCode code = DindelUtil::runDindelPair(refStr,
                                                      varStr,
                                                      m_parameters,
                                                      m_baseVCFFile,
                                                      m_variantVCFFile);

    m_returnCodes[code] += 1;
    
}

std::string VCFTester::applyVariant(const std::string& in, int pos,
                                    const std::string& ref, const std::string& alt)
{
    std::string out = in;
    assert(pos > 0 && pos < (int)out.length());

    // Ensure that the reference string at the variant matches the expected
    assert(out.substr(pos, ref.length()) == ref);
    out.replace(pos, ref.length(), alt);
    return out;
}
