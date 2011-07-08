//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFUtil - collection of utility functions
// for producing VCF records.
//
#include <iostream>
#include <stdlib.h>
#include <iterator>
#include "VCFUtil.h"
#include "StdAlnTools.h"
#include "MultiAlignment.h"

// 
// VCFRecord
//
void VCFRecord::addComment(const std::string& key, const std::string& value)
{
    std::string out = key;
    out.append(1, '=');
    out.append(value);


    comments.push_back(out);
}

void VCFRecord::printVerbose() const
{
    std::cout << "RefName: " << refName << "\n";
    std::cout << "Position: " << refPosition << "\n";
    std::cout << "RefStr: " << refStr << "\n";
    std::cout << "VarStr: " << varStr << "\n";
}

// Write the VCF record to the stream
std::ostream& operator<<(std::ostream& o, VCFRecord& record)
{
    o << record.refName << "\t";
    o << record.refPosition << "\t",
    o << ".\t"; // ID not supported
    o << record.refStr << "\t";
    o << record.varStr << "\t";
    o << ".\t"; // QUAL not supported
    o << "PASS\t"; // everything passes
    std::copy( record.comments.begin(), record.comments.end(), 
               std::ostream_iterator<std::string>(o, ";")); // comments 
    return o;
}


//
VCFClassification VCFRecord::classify() const
{
    if(refStr.size() == varStr.size())
        return VCF_SUB;
    
    if(refStr.size() > varStr.size())
    {
        // If a deletion, varStr should be a prefix of refStr
        if(refStr.find(varStr) == 0)
            return VCF_DEL;
        else
            return VCF_COMPLEX;
    }

    if(refStr.size() < varStr.size())
    {
        // If an insertion, refStr should be a prefix of refStr
        if(varStr.find(refStr) == 0)
            return VCF_INS;
        else
            return VCF_COMPLEX;
    }
    return VCF_COMPLEX;
}

//
// VCFUtil
//

// See header for description
VCFReturnCode VCFUtil::generateVCFFromCancerVariant(const std::string& ref,
                                                    const std::string& base,
                                                    const std::string& variant,
                                                    const std::string& refName,
                                                    size_t refPosition,
                                                    const std::string& varName,
                                                    int minExactMatch,
                                                    VCFVector& outRecords)
{
    VCFVector tempRecords;

    // Construct a multiple alignment from pairwise alignments
    // of the base and variant strings to the reference

    // We start by constructing a global alignment between the sequences
    // and representing the alignment as a Cigar string
    std::string cigarRefBase = StdAlnTools::globalAlignmentCigar(base, ref);
    std::string cigarRefVariant = StdAlnTools::globalAlignmentCigar(variant, ref);
    
    // Initialize the multiple alignment data
    
    MAlignData baseMA;
    baseMA.position = 0;
    baseMA.str = base;
    baseMA.expandedCigar = StdAlnTools::expandCigar(cigarRefBase);
    baseMA.name = "base";

    MAlignData varMA;
    varMA.position = 0;
    varMA.str = variant;
    varMA.expandedCigar = StdAlnTools::expandCigar(cigarRefVariant);
    varMA.name = "variant";

    MAlignDataVector maVector;
    maVector.push_back(baseMA);
    maVector.push_back(varMA);
    
    MultiAlignment ma(ref, maVector);

    std::cout << "\nMultipleAlignment:\n";
    ma.print();

    // Get the row index of the variant and base strings
    size_t baseRow = ma.getIdxByName("base");
    size_t varRow = ma.getIdxByName("variant");
    size_t refRow = ma.getRootIdx();

    size_t numCols = ma.getNumColumns();

    // Iterate over every column of the multiple alignment and detect variants
    int leftExactMatch = 0;
    int rightExactMatch = 0;

    int eventStart = -1;
    int eventEnd = -1;
    bool isIndel = false;

    bool inLeftExact = true;

    for(size_t i = 0; i < numCols; ++i)
    {
        char baseSymbol = ma.getSymbol(baseRow, i);
        char varSymbol = ma.getSymbol(varRow, i);
        bool isVariant = (varSymbol != baseSymbol);

        // Update the counter of the number of exact matches for the leftmost and rightmost bases
        if(isVariant)
        {
            inLeftExact = false; // stop the count of leftmost exact matches
            rightExactMatch  = 0; // reset the rightmost count
        }
        else // this is a good match
        {
            rightExactMatch += 1;
            if(inLeftExact)
                leftExactMatch += 1;
        }
        
        // Process variants indicated by this portion of the alignment
        if(isVariant)
        {
            if(eventStart == -1)
            {
                // Not currently in a variant event, start a new one
                eventStart = i;
                eventEnd = i;
            }
            else
            {
                // Extend the event
                assert(eventEnd != -1);
                eventEnd += 1;
            }

            isIndel = isIndel || varSymbol == '-' || baseSymbol == '-';
        }
        else
        {
            // Check if this is the end of a variant
            if(eventStart != -1)
            {
                // Extract the substrings of the reference and variant.
                // If this is an indel event, we start the base one base
                // upstream of the variant
                if(isIndel)
                {
                    if(eventStart == 0)
                    {
                        std::cerr << "VCFUtil error: indel started at position 0\n";
                        exit(EXIT_FAILURE);
                    }   
                    eventStart -= 1;
                }

                size_t eventLength = eventEnd - eventStart + 1;
                std::string refString = StdAlnTools::unpad(ma.getPaddedSubstr(refRow, eventStart, eventLength));
                std::string varString = StdAlnTools::unpad(ma.getPaddedSubstr(varRow, eventStart, eventLength));
                
                // Generate VCF Record
                VCFRecord record;
                record.refName = refName;
                record.refPosition = refPosition + eventStart + 1; // VCF uses 1-based coordinates
                WARN_ONCE("DONT USE COL - USE REFERENCE BASE IDX");
                printf("Ref start: %d col: %d eventStart: %d vcf pos: %d\n", (int)refPosition, (int)i, (int)eventStart, (int)record.refPosition);
                record.refStr = refString;
                record.varStr = varString;
                record.addComment("id", varName);
                tempRecords.push_back(record);

                // Reset state
                eventStart = -1;
                eventEnd = -1;
                isIndel = false;
            }
        }
    }

    // Perform the start/end exact match sanity check
    if(leftExactMatch >= minExactMatch && rightExactMatch >= minExactMatch)
    {
        outRecords.insert(outRecords.end(), tempRecords.begin(), tempRecords.end());
        return VCF_OK;
    }
    else
    {
        return VCF_EXACT_MATCH_FAILED;
    }
}

// Write a VCF header to the given file handle
void VCFUtil::writeHeader(std::ostream* pWriter, const std::string& program, const std::string& inputBAM, const std::string& reference)
{
    // Get and format the current date
    time_t rawtime;
    struct tm * timeinfo;
    static int MAX_DATE_CHARS = 80;
    char dateBuffer [MAX_DATE_CHARS];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    int copied = strftime(dateBuffer,MAX_DATE_CHARS,"%Y%m%d",timeinfo);
    assert(copied < MAX_DATE_CHARS);

    *pWriter << "##fileformat=VCFv4.0\n";
    *pWriter << "##fileDate=" << dateBuffer << "\n";
    *pWriter << "##command=" << program << "\n";
    *pWriter << "##inputBAM=" <<  inputBAM << "\n";
    *pWriter << "##reference=" << reference << "\n";
    *pWriter << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}
