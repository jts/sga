//-----------------------------------------------
// Copyright 201- Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
//        and Kees Albers (caa@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DindelRealignWindow - Infer haplotypes using the FM-index
//
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>
#include <iomanip>
#include "StdAlnTools.h"
#include "MultiAlignment.h"
#include "DindelRealignWindow.h"
#include "BWT.h"
#include "DindelHMM.h"
#include <cmath>

const int DINDEL_DEBUG=0;
const int QUIET=1;
const int ADDSNPS=1;
const int ALWAYS_REALIGN=1;
const int DEBUG_CALLINDEL=0;
const int REPOSITION_INDEL_WINDOW=1000;
//#define OVERLAPPER // build overlapper

#include <iostream>


long long int combinations(const int n, const int k)
{
    // If n equals with k or k equals to 0, then the result is always 1
    if (n == k || k == 0)
        return 1;
    // If k is 1, then the result is always n
    if (k == 1)
        return n;

    // n < k, can't calculate, return -1
    if (n < k)
        return -1;
    else
    {
        int i;
        long long int result, fact_n, fact_k, fact_n_sub_k;

        // Calculate n!
        fact_n = 1;
        for (i = 2; i <= n; i++)
            fact_n *= i;

        // Calculate k!
        fact_k = 1;
        for (i = 2; i <= k; i++)
            fact_k *= i;

        // Calculate (n - k)!
        fact_n_sub_k = 1;
        i = 2;
        while (i <= (n - k))
        {
            fact_n_sub_k *= i;
            i++;
        }

        result = fact_n / (fact_k * fact_n_sub_k);
        return result;
    }
}


void parseRegionString(const std::string & region, std::string & chrom, int & start, int & end)
{
        std::string filtered;
        for(size_t x=0;x<region.size();x++) {
                char c=region[x];
                if (c=='-' || c==':') filtered+=' ';
                else filtered+=c;
        }
        std::istringstream is(filtered);
        is >> chrom;
        std::string e; is >> e;
        if (!from_string(start,e,std::dec)) throw std::string("Cannot parse region for start!");
        is >> e;
        if (!from_string(end,e,std::dec)) throw std::string("Cannot parse region end!");

        if (end<start)
        {
            std::string message("Invalid region: start>end : ");
            message+=region;
            throw message;
        }
}

std::vector<std::string> SplitString(const std::string & str, char sep)
{
    std::string elem;
    std::vector<std::string> split;

    for (size_t x=0;x<str.size();x++)
    {
        if (str[x]==sep)
        {
            if (!elem.empty()) split.push_back(elem);
            elem = "";
        } else elem += str[x];

    }
    if (!elem.empty()) split.push_back(elem);
    return split;
}

/*

 * DINDELBAM

 */








/*
 *
 * DINDELVARIANT
 *
 *
 */

DindelRead::DindelRead(const SeqItem & seqItem, const SampleName & sampleName, double mappingQual, int fixedBaseQual, bool isForward) : m_seqItem(seqItem)
{
    m_mappingQual = mappingQual;
    m_fixedBaseQual = fixedBaseQual;
    m_rcRead = false;
    m_sampleName = sampleName;
    m_setupHash = false;
    m_isForward = isForward;
}

void DindelRead::getLogProbCorrectError(std::vector<double>& lpCorrect, std::vector<double>& lpError) const
{
    lpError=std::vector<double>(this->length());
    lpCorrect=std::vector<double>(this->length());
    for (int b=0;b<this->length();b++)
    {
        double pr=1.0 - exp ( (-2.3026/10.0)*double(this->getQual(b) ));
        double eq=log(.25+.75*pr);
        double uq=log(.75+1e-10-.75*pr);
        lpError[b] = uq;
        lpCorrect[b] = eq;
    }
}


void DindelRead::setupHash()
{
    if (!m_rcRead)
        DindelSequenceHash::getKeys(m_hashKeys, this->getSequence());
    else
    {
        std::string seq = this->getSequence();
        std::string rcSeq=seq;
        size_t rlen=rcSeq.size();
        for (size_t x=0;x<rcSeq.size();x++) rcSeq[rlen-x-1]=_complement(seq[x]);
        DindelSequenceHash::getKeys(m_hashKeys, rcSeq);
    }
    m_setupHash = true;
}

void DindelSequenceHash::getKeys(std::vector<unsigned int>& keys, const std::string& sequence)
{

    if (sequence.size()>=DINDEL_HASH_SIZE)
    {
        keys.clear();
        keys.reserve(sequence.size()-DINDEL_HASH_SIZE);
        unsigned int key=DindelSequenceHash::convert(sequence, 0);
    keys.push_back( key );
        for (size_t x=1;x<sequence.size()-DINDEL_HASH_SIZE;x++)
        {
        key = DindelSequenceHash::pushBack(key, sequence[x+DINDEL_HASH_SIZE]);
            keys.push_back( key );

        }
    }    
}

const std::vector<unsigned int> & DindelRead::getHashKeys()
{
    if (!m_setupHash) setupHash();
    return m_hashKeys;
}

VariantPriors::VariantPriors()
{
    m_probSNP = 0.001;
    m_probINDEL = 0.0001;
    m_probMNP = 0.00001;
}


double VariantPriors::getDefaultProbVariant(const std::string& type) const
{
    if (type == "SNP") return m_probSNP;
    else if (type == "INDEL") return m_probINDEL;
    else if (type == "MNP") return m_probMNP;
    else assert(1==0); // ("Unknown variant type");
}


DindelVariant::DindelVariant()
{
    m_pos = -1;
    m_type = "uninitialized";
}

DindelVariant::DindelVariant(const std::string & chrom, const std::string & ref, const std::string & alt, int pos)
{
    if (ref.empty() || alt.empty() || chrom.empty()) throw std::string("DindelVariant::zero_length_chrom_ref_alt_string");
    if (pos<0) throw std::string("DindelVariant::negative_variant_position");

    m_chrom = chrom;
    m_ref = ref;
    m_alt = alt;
    m_pos = pos;
    
    m_leftUniquePos = -1;
    m_rightUniquePos = -1;
    m_idx = -1;
    m_priorProb=-1;
    
    setup();
}


bool DindelVariant::variantFromWindowString(const std::string & varstring, DindelVariant & variant)
{
    // split string

    std::vector<std::string> vs;
    size_t lastPos = 0;
    for (size_t x=0;x<varstring.size();x++)
    {
        if (varstring[x]==',')
        {
            if (x>0)
            {
                vs.push_back(varstring.substr(lastPos, x-lastPos));
                lastPos=x+1;
            }
        }
    }
    if (lastPos<varstring.size())
    {
        vs.push_back(varstring.substr(lastPos, varstring.size()-lastPos));
    }

    if (vs.size()<3) return false;

    int pos;
    std::string chrom, ref, alt;
   
    chrom = vs[0]; if (chrom.empty()) return false;
    if (!from_string<int>(pos, vs[1], std::dec)) return false;
    if (pos<0) return false;

    ref = vs[2]; if (ref.empty()) return false;
    alt = vs[3]; if (alt.empty()) return false;

    for (size_t x=0;x<ref.size();x++) if (ref[x]!='A' && ref[x]!='T' && ref[x] != 'C' && ref[x]!='G' && ref[x] != 'N') return false;
    for (size_t x=0;x<alt.size();x++) if (alt[x]!='A' && alt[x]!='T' && alt[x] != 'C' && alt[x]!='G' && alt[x] != 'N') return false;

    variant = DindelVariant(chrom, ref, alt, pos);

    return true;
}


void DindelVariant::checkSequence(const std::string & seq)
{
    for (size_t x=0;x<seq.size();x++) if (seq[x]!='A' && seq[x] !='C' && seq[x] != 'G' && seq[x] != 'T' && seq[x] != 'N') throw std::string("Variant: sequence "+seq+" not allowed.");

}

void DindelVariant::setup()
{
    // check if there are no commas in alt or ref

    checkSequence(m_alt);
    checkSequence(m_ref);

    m_dlen = m_alt.length()-m_ref.length();

    if (m_dlen==0 && m_alt.size()==1) m_type = "SNP";
    else if (m_dlen==0 && m_alt.size()>1) m_type = "MNP";
    else if (m_dlen!=0) m_type = "INDEL";

    std::stringstream posStr;
    posStr << m_pos;

    m_id = m_chrom+"_"+posStr.str()+"_"+m_ref+"_"+m_alt;
}

int DindelVariant::getHaplotypeLeftUnique() const
{
    if (m_leftUniquePos==-1) throw std::string("DindelVariant::call_determineLeftRightUnique_first");
    return m_leftUniquePos;
}

int DindelVariant::getHaplotypeRightUnique() const
{
    if (m_rightUniquePos==-1) throw std::string("DindelVariant::call_determineLeftRightUnique_first");
    return m_rightUniquePos;
}

void DindelVariant::write(std::ostream & out) const
{
    out << m_id << " " << m_chrom << " " << m_ref << " " << m_alt << " type: " << m_type << " priorProb: " << m_priorProb << std::endl;
}

/*
 *
 * DINDELHAPLOTYPE
 *
 *
 */

DindelHaplotype::DindelHaplotype(const std::string & refSeq, int refSeqStart, bool isReference = false)
{
    assert(isReference == true);
    m_seq = refSeq;
    m_hapRefStart = refSeqStart;
    m_isReference = isReference; // it is only a reference haplotype if it is explicitly indicated
    initHaplotype();
}

DindelHaplotype::DindelHaplotype(const std::string & refName, const std::string & refSeq, int refSeqStart, const MultiAlignment & ma, size_t varRow, size_t refRow)
{
    size_t numCols = ma.getNumColumns();
    if (DINDEL_DEBUG) std::cout << "DindelHaplotype::DindelHaplotype numCols: " << numCols << std::endl;

    

    this->m_seq = StdAlnTools::unpad(ma.getPaddedSubstr(varRow,0,numCols));
    if (DINDEL_DEBUG)
    {
        ma.print();
        std::cout << "REFSEQ: " << refSeq << "\n";
        std::cout << "UNPAD:  " << StdAlnTools::unpad(ma.getPaddedSubstr(refRow,0,numCols)) << "\n";
    }
    assert(refSeq == StdAlnTools::unpad(ma.getPaddedSubstr(refRow,0,numCols)));

    

    m_refPos = std::vector<int>(m_seq.size(), 0);
    determineHomopolymerLengths();
    m_sequenceHash=DindelSequenceHash(m_seq);

    this->m_isReference = false;


    // set m_refPos
    int hidx=-1, ridx=-1, hDelStart=-1;
    size_t numLeftOverhang = 0, numRightOverhang = 0;
    bool inRefDeletion = false;
    bool leftOverhang = true;
     
    /*
      note that the following case still needs to be dealt with if left or right overhang is allowed
    0       TGCTATTCTCTCCAACAAGACCGTTGAAC----------A        REF.chr10
    1       TGCTATTCTCTCCAACAAGACCGTTGAACAATTGGGGCAA        haplotype-1
    2       TGCTATTCTCTCCAACAAGACCGTTGAACAATTGGGGCAA        haplotype-2
    */

    for(size_t i = 0; i < numCols; ++i)
    {
        char rs = ma.getSymbol(refRow,i);
        char vs = ma.getSymbol(varRow,i);
        
        if (rs != '-') ridx++;
        if (vs != '-') hidx++;

        if (leftOverhang && rs != '-')
        {
            // end of left overhang of haplotype with reference
            numLeftOverhang = i;
            for (size_t j=0;j<i;j++) m_refPos[j]= LEFTOVERHANG;
            leftOverhang = false;
        }
        if (!leftOverhang)
        {
            if (rs == '-' && vs != '-') m_refPos[hidx] = INSERTION;
            if (rs != '-' && vs != '-')
            {
                if(rs!=vs) m_refPos[hidx] = SNP; else m_refPos[hidx] = refSeqStart+ma.getBaseIdx(refRow, i);

            }
            if (rs == '-' && !inRefDeletion)
            {
                inRefDeletion = true;
                hDelStart = hidx;
            }
            if (rs != '-') inRefDeletion = false;
        }
    }

    if (inRefDeletion)
    {
        // ref deletion extends to end of alignment, right overhang.
        assert(hDelStart>0);
        assert(m_refPos[hDelStart-1] >=0);
        for(int h = hDelStart; h < int(m_refPos.size()); h++)
        {
            m_refPos[h] = RIGHTOVERHANG;
            numRightOverhang++;
        }
    }

    if(numLeftOverhang > 0 || numRightOverhang > 0)
        throw std::string("DindelHaplotype::DindelHaplotype No overhanging bases");
    //assert(numLeftOverhang == 0 && numRightOverhang == 0);

    if(DINDEL_DEBUG)
    {
        std::cout << "m_refPos: ";
        for (size_t i = 0; i < m_refPos.size(); i++) std::cout << "[" << i << " " << m_refPos[i] << "]";
        std::cout << "\n";
    }

    int verbose = 0;

    int leftExactMatch = 0;
    int rightExactMatch = 0;

    int eventStart = -1;
    int eventEnd = -1;
    bool isIndel = false;
    bool isComplex = false;

    bool inLeftExact = true;

    for(size_t i = 0; i < numCols; ++i)
    {
        char refSymbol = ma.getSymbol(refRow, i);
        char varSymbol = ma.getSymbol(varRow, i);

        // JS 25/10/11 Keep considering the event to be a variant if the var or ref
        // sequence has a gap here
        bool isVariant = (varSymbol != refSymbol || varSymbol == '-' || refSymbol == '-');

        if (DINDEL_DEBUG) std::cout << " ** i: " << i << " isVariant: " << isVariant << " refSymbol: " << refSymbol << " varSymbol: " << varSymbol << std::endl;


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
                // NOTE that a SNP next to an indel will be annotated as a single event.
                isComplex = (!isIndel && (varSymbol == '-' || refSymbol == '-')) || (isIndel && varSymbol !='-' && refSymbol != '-');
                assert(eventEnd != -1);
                eventEnd += 1;
            }

            // event will be an indel even if it contains single nucleotide mismatches to the reference

            isIndel = isIndel || varSymbol == '-' || refSymbol == '-';

        }
        else
        {
            // Check if this is the end of a variant
            if(eventStart != -1)
            {
                if(verbose > 0)
                {
                    std::cout << "\n\nProcessing variant\n";
                    ma.print();
                }

                // Extract the substrings of the three sequences describing
                // this event.
                // If the event is an indel, we move the start point
                // of the event backwards until all 3 columns are not padding characters.
                // This is to avoid the following situation:
                // AATATTT--CGT ref
                // AATATTTT-CGT base
                // AATATTTTTCGT var
                // or
                // TTAGGTTTTTTT ref
                // TTAGG-TTTTTT base
                // TTAGG--TTTTT var
                //
                // In the first case, it is a 1 base insertion wrt base but a 2 base insertion wrt to the ref.
                // So we need to move the event to the first column with all Ts to properly capture what is happening.
                // In the second case it is a 1 base deletion wrt base and a 2 base deletion wrt to the ref.
                // Again, we need to move back to the last full column.

                bool okEvent = true;
                if(isIndel && (eventStart == 0 || eventEnd == (int)numCols - 1) ) okEvent = false;

                if(okEvent)
                {
                    if(isIndel)
                    {
                        // It is impossible to extract a reference string if the indel is at the beginning
                        // or end of th multiple alignment so we just fail out here

                        while(eventStart >= 0)
                        {
                            char rc = ma.getSymbol(refRow, eventStart);
                            char vc = ma.getSymbol(varRow, eventStart);
                            if(rc == '-' || vc == '-')
                                eventStart -= 1;
                            else
                                break;
                        }

                        assert(eventStart >= 0);
                    }

                    size_t eventLength = eventEnd - eventStart + 1;
                    std::string refStringPadded = ma.getPaddedSubstr(refRow, eventStart, eventLength);
                    std::string varStringPadded = ma.getPaddedSubstr(varRow, eventStart, eventLength);
                    std::string refString = StdAlnTools::unpad(refStringPadded);
                    std::string varString = StdAlnTools::unpad(varStringPadded);
                    int hap_start = int(ma.getBaseIdx(varRow, eventStart));
                    int hap_end = int(ma.getBaseIdx(varRow,eventStart+eventLength));
                    if(verbose > 0)
                    {
                        printf("Ref start: %d col: %d eventStart: %d eventEnd: %d\n", (int)refSeqStart, (int)i, eventStart, eventEnd);
                        std::cout << "RefString " << refString << "\n";

                        std::cout << "varString " << varString << "\n";
                    }
                    assert(!refString.empty());
                    //assert(!baseString.empty());
                    assert(!varString.empty());

                    // minPos and maxPos are coordinates of window in MultiAlignment where variant can be ambiguously positioned
                    int minPos = eventStart, maxPos = eventStart;

                    if(isIndel)
                    {
                        // If it is a clean indel event, see how far it can be moved to the left or right
                        bool isDeletion = refString.size() > varString.size();
                        std::string refPadded = ma.getPaddedSubstr(refRow,0,numCols);
                        int indelLength = int(eventLength) - 1;
                        assert(indelLength>=1);

                        // ---AT-GAATCG-T-- REF  There might be multiple events
                        // ATCATCGA--CGATCG ALT

                        // ignore overlaps of the haplotype with the reference at the ends

                        minPos = refPadded.size()-1;
                        maxPos = 0;

                        // create candidate haplotype with just this variant removed
                        

                        std::string changeHaplotype = ma.getPaddedSubstr(varRow,0,numCols);
                        changeHaplotype.replace(eventStart, eventLength, ma.getPaddedSubstr(refRow, eventStart, eventLength));
                        std::string indelSeq = isDeletion ? refStringPadded.substr(1,indelLength) : varStringPadded.substr(1,indelLength) ;

                        if (DINDEL_DEBUG) std::cout << "CHANGEHAPLOTYPE: " << changeHaplotype << std::endl;
                        

                        for(int j=0;j<int(numCols)-indelLength;j++) if (ma.getSymbol(refRow, j)!='-')
                        {
                            std::string alt(changeHaplotype);

                            int k=0,numModified=0;
                            if (isDeletion)
                            {
                                while (numModified<indelLength && j+k<int(numCols))
                                {
                                    if(alt[j+k]!='-')
                                    {
                                        alt[j+k]='-';
                                        numModified++;
                                    }
                                    k++;
                                }

                            } else
                            {
                                alt.insert(size_t(j),indelSeq);
                                numModified=indelLength;
                            }

                            std::string altUnpadded = StdAlnTools::unpad(alt);
                            if(numModified == indelLength && altUnpadded == m_seq)
                            {
                                if (DINDEL_DEBUG) std::cout << "CHANGEHAPLOTYPE equal " << j << std::endl;
                                if(j<minPos) minPos = j;
                                if(j>maxPos) maxPos = j;
                            }
                        }


                    }

                    // Get the base position in the reference string. This is not necessarily the
                    // same as the eventStart column as the reference may be padded
                    // int refBaseOffset = ma.getBaseIdx(refRow, eventStart);
                    // get positions in haplotype where indel variant may be ambiguously positioned

                    int refBaseOffsetMinPos = ma.getBaseIdx(refRow, minPos);
                    

                    while (minPos>0 && ma.getSymbol(varRow, minPos)== '-' ) minPos--;
                    while (maxPos<int(numCols)-1 && ma.getSymbol(varRow, maxPos)== '-' ) maxPos++;

                    int varBaseOffsetMinPos = ma.getBaseIdx(varRow, minPos);
                    int varBaseOffsetMaxPos = ma.getBaseIdx(varRow, maxPos);

                    // Use leftmost position for the variant (which can only be change for an insertion or deletion)
                    DindelVariant var(refName, refString, varString, refBaseOffsetMinPos+refSeqStart);
                    var.setPriorProb(0.001); //FIXME
                    if (DINDEL_DEBUG) std::cout << "VARIANT: " << refName << " " << refString << "/" << varString << " pos: " << refBaseOffsetMinPos+refSeqStart << " varBaseOffsetMinPos: " << varBaseOffsetMinPos << " varBaseOffsetMaxPos: " << varBaseOffsetMaxPos << std::endl;
                    var.setHaplotypeUnique(varBaseOffsetMinPos, varBaseOffsetMaxPos);
                    m_variants.push_back(var);

                    std::pair < HashMap<std::string, std::pair<int, int> >::iterator, bool> ins_pair =  m_variant_to_pos.insert( HashMap<std::string, std::pair<int, int> >::value_type ( var.getID(), std::pair<int,int>(hap_start, hap_end)));
                    assert (ins_pair.second == true);
                }
                // Reset state
                eventStart = -1;
                eventEnd = -1;
                isIndel = false;
                isComplex = false;
            }
        }



    }
    



}

DindelHaplotype::DindelHaplotype(const DindelHaplotype & haplotype, int copyOptions)
{
   if (copyOptions == 0) 
   {
       // copy variants and m_refPos from haplotype
       
       m_seq = haplotype.m_seq;
       m_refPos = haplotype.m_refPos;
       m_hplen = haplotype.m_hplen;
       m_variants = haplotype.m_variants;
       m_variant_to_pos = haplotype.m_variant_to_pos;
       m_hapRefStart = haplotype.m_hapRefStart;
       m_isReference = haplotype.m_isReference; 
       // DO NOT CALL initHaploype()
   }
}

void DindelHaplotype::initHaplotype()
{
    m_refPos = std::vector<int>(m_seq.size(), 0);

    // set reference positions
    for (size_t r=0;r<m_seq.size();r++)
    {
        m_refPos[r] = m_hapRefStart+int(r);
    }

    determineHomopolymerLengths();
    m_sequenceHash=DindelSequenceHash(m_seq);
}

void DindelHaplotype::updateHaplotype()
{

    determineHomopolymerLengths();

    // check if the first haplotype base is still aligned to the reference sequence
    m_hapRefStart = m_refPos[0];
    if (m_hapRefStart<0) throw std::string("DindelHaplotype::updateHaplotype::first_haplotype_base_not_on_reference");

    m_sequenceHash=DindelSequenceHash(m_seq);
}


void DindelHaplotype::determineHomopolymerLengths()
{
   
    // determine homopolymer length in sequence

    m_hplen = std::vector<int>(m_seq.size(), 1);

    std::vector<int> hp_forward = std::vector<int>(m_seq.size(), 0);
    std::vector<int> hp_reverse = std::vector<int>(m_seq.size(), 0);

    for (size_t r=1;r<m_seq.size();r++)
    {
        int add=0;
        if (m_seq[r]==m_seq[r-1]) add=hp_forward[r-1]+1;
        hp_forward[r] = add;
    }

    for (int r=int(m_seq.size())-2;r>=0;r--)
    {
        int add=0;
        if (m_seq[r]==m_seq[r+1]) add=hp_reverse[r+1]+1;
        hp_reverse[r] = add;
    }

    for (size_t r=0;r<m_seq.size();r++) m_hplen[r] = hp_reverse[r]+hp_forward[r]+1;
}

void DindelHaplotype::write(std::ostream & out) const
{
    out << "VARIANTS: "; for (size_t x=0;x<m_variants.size();x++) out << " " << m_variants[x].getID(); out << std::endl;
    for (size_t x=0;x<m_seq.size(); x++) out << m_seq[x]; out << std::endl;
    for (size_t x=0;x<m_refPos.size(); x++) out << " " << m_refPos[x]-m_hapRefStart; out << std::endl;
}

bool DindelHaplotype::addVariant(const DindelVariant& var)
{
    if (DINDEL_DEBUG) std::cerr << "DindelHaplotype::addVariant START" << std::endl;
    // returns true if the variant was added
    int relpos = getHapBase(var.getPos());

    if (relpos<10 || relpos>this->length()-10) return false; // don't add it if it would change the first base of the haplotype

    // check if all ref bases in haplotype that var wants to change are still marked as reference
    bool refChanged = false;
    for (int p=relpos;p<relpos+int(var.getAlt().size());p++)
    {
        if (m_refPos[p]<0)
        {
            refChanged = true;
            break;
        }
    }

    // check whether the variant being added is a SNP and is not too close to an indel

    if (var.getType()=="SNP" || var.getType()=="MNP")
    {
        int l=relpos-5; if (l<0) l=0;
        int r=relpos+5; if (r>this->length()) r=this->length();
        for (int b=l;b<r;b++)
        {
            if (m_refPos[b]==DELETION_NOVEL_SEQUENCE || m_refPos[b]==INSERTION)
            {
                refChanged = true; // don't add the SNP, it is too close to an indel.
                break;
            }
        }

    }

    if (!refChanged) {
        // the reference is still there, let's change the haplotype
        // update m_seq and m_refPos

        const std::string & ref = var.getRef();
        const std::string & alt = var.getAlt();
        int refLen = (int) ref.length();
        int altLen = (int) alt.length();
        int dlen = var.getLengthDifference();
        // determine number of bases that are identical starting from the left and the rights

        int leftIdentical = 0;
        int rightIdentical = 0;

        for (int x=0;x<altLen;x++)
        {
            if (alt[x]==ref[x]) leftIdentical++; else break;
        }

        for (int x=0;x<altLen;x++)
        {
            if (alt[altLen-x-1]==ref[refLen-x-1]) rightIdentical++; else break;
        }
        

        if (leftIdentical+rightIdentical>altLen)
        {
            leftIdentical=altLen-rightIdentical; // this will shift the first variant base to the left
        }

        // change the haplotype
    assert(m_seq.substr(relpos, refLen)==ref);
        m_seq.replace(relpos, refLen, alt);

        if (DINDEL_DEBUG) std::cerr << "HAPADDVARIANT " << var.getPos() << " " << ref << " " << alt << "relpos: " << relpos << " leftIdentical: " << leftIdentical << " rightIdentical: " << rightIdentical << " refLen: " << refLen << " altLen: " << altLen << std::endl;

        // update m_refPos
        int hap_start, hap_end;
        if (dlen<0)
        {
            // some bases from the reference will be deleted
            std::vector<int> tmp(m_seq.size(),ADDVARIANT_DEBUG);
            for (int r=0;r<relpos+leftIdentical;r++) tmp[r] = m_refPos[r];

            // these are base that were not in the reference sequence
            for (int r=relpos+leftIdentical;r<relpos+altLen-rightIdentical;r++) tmp[r] = DELETION_NOVEL_SEQUENCE;

            // takes care of the bit on the right that is identical
            for (int r=relpos+altLen-rightIdentical;r<relpos+altLen;r++) tmp[r] = m_refPos[r-dlen];

            // this takes care of the bit after the change
            for (int r=relpos+altLen;r<int(m_seq.size());r++) tmp[r] = m_refPos[r-dlen];

            m_refPos.swap(tmp);

            hap_start = relpos+leftIdentical;
            hap_end = relpos+altLen-rightIdentical;

        } else if (dlen>0)
        {
            // some bases will be added

            std::vector<int> tmp(m_seq.size(),ADDVARIANT_DEBUG);
            for (int r=0;r<relpos+leftIdentical;r++) tmp[r] = m_refPos[r];

             // these are base that were not in the reference sequence
            for (int r=relpos+leftIdentical;r<relpos+altLen-rightIdentical;r++) tmp[r] = INSERTION;
            // takes care of the bit on the right that is identical
            for (int r=relpos+altLen-rightIdentical;r<relpos+altLen;r++) tmp[r] = m_refPos[r-dlen];


            // this takes care of the bit after the change
            for (int r=relpos+altLen;r<int(m_seq.size());r++) tmp[r] = m_refPos[r-dlen];

            m_refPos.swap(tmp);

            hap_start = relpos+leftIdentical;
            hap_end = relpos+altLen-rightIdentical;
        if (DINDEL_DEBUG) {
            std::cerr << " hap_start: "<< hap_start << std::endl;
            std::cerr << " hap_end: " << hap_end << std::endl;
        }


        } else {
            // length of sequence is unchanged. Might be multinucleotide SNP
        if (DINDEL_DEBUG) std::cerr << "ADDING SNP" << std::endl;
            for (int r=relpos+leftIdentical;r<=relpos+rightIdentical;r++) m_refPos[r] = MULTINUCLEOTIDE_RUN;
            hap_start = relpos+leftIdentical;
            hap_end = relpos+rightIdentical;
        }

        // check if m_refPos and m_seq are still consistent
        if (m_seq.size() != m_refPos.size()) 
    {
        std::cerr << "HAPADDVARIANT " << var.getPos() << " " << ref << " " << alt << " relpos: " << relpos << " leftIdentical: " << leftIdentical << " rightIdentical: " << rightIdentical << " refLen: " << refLen << " altLen: " << altLen << " haplotype length: " << this->length() << std::endl;   
        this->write(std::cerr);
        throw std::string("DindelHaplotype::addVariant::m_seq_and_m_refPos_inconsistent");
    }
        for (size_t x=0;x<m_refPos.size();x++) if (m_refPos[x]==ADDVARIANT_DEBUG) throw std::string("DindelHaplotype::addVariant::m_refPos_error");

        if (dlen>0)
        {
            bool hasIns = false;
            for (size_t x=0;x<m_refPos.size();x++) if (m_refPos[x]==INSERTION) hasIns = true;
            if (!hasIns) throw std::string("DindelHaplotype::addVariant::m_seq_and_m_refPos_inconsistent_inserror");
        } else if (dlen<0)
        {
            bool hasDel = false;
            for (size_t x=1;x<m_refPos.size();x++) if (m_refPos[x]-m_refPos[x-1]>1) hasDel = true;
            if (!hasDel) throw std::string("DindelHaplotype::addVariant::m_seq_and_m_refPos_inconsistent_delerror");
        }

    m_variants.push_back(var);

    std::pair < HashMap<std::string, std::pair<int, int> >::iterator, bool> ins_pair =  m_variant_to_pos.insert( HashMap<std::string, std::pair<int, int> >::value_type ( var.getID(), std::pair<int,int>(hap_start, hap_end)));

        // set coordinates of variant in haplotype
    //assert(this->m_variant_to_pos.find(var.getID())!=this->m_variant_to_pos.end());

    m_isReference=false;
        updateHaplotype();
        if (DINDEL_DEBUG) std::cerr << "DindelHaplotype::addVariant ADDED END" << std::endl;

        return true;
    } else
    {
    if (DINDEL_DEBUG) std::cerr << "DindelHaplotype::addVariant **NOT** ADDED END" << std::endl;    
        return false;
    }

}

int DindelHaplotype::getClosestDistance(const DindelVariant& variant, int hapPos1, int hapPos2) const
{
    HashMap<std::string, std::pair<int, int> >::const_iterator it = this->m_variant_to_pos.find(variant.getID());
    if (it == this->m_variant_to_pos.end())
    {
    std::cerr << " query variant: " << variant.getID() << std::endl;
    std::cerr << " haplotype variants: "; for (it = m_variant_to_pos.begin();it != m_variant_to_pos.end();it++) std::cerr << " " << it->first; std::cerr << std::endl;
    std::cerr << " h2 variants:        "; for (size_t x=0;x<m_variants.size();x++) std::cerr << " " << m_variants[x].getID(); std::cerr << std::endl;
        return -1;
    } else {
        int s = it->second.first;
        int e = it->second.second;
    
    // check if there is any overlap at all.
    if (hapPos1>hapPos2) 
    {
        int t = hapPos1;
        hapPos1 = hapPos2;
        hapPos2 = t;
    }
    
    if (hapPos1<=s && hapPos2>=e) 
    {
        std::set<int> dists;
        /*
        dists.insert(abs(hapPos1-s));
        dists.insert(abs(hapPos2-s));
        dists.insert(abs(hapPos1-e));
        dists.insert(abs(hapPos2-e));
        */
        // make sure we only count extensions from the left and right
        int d1=hapPos2-e;
        if (d1<0) d1=0;
        int d2=s-hapPos1;
        if (d2<0) d2=0;
        dists.insert( d1 );
        dists.insert( d2 );

        return *dists.begin();
    } else 
    {
        // no overlap
        return -1;
    }
    }
}

int DindelHaplotype::getHapBase(int refPos) const
{
    for (size_t x=0;x<this->m_refPos.size();x++) if (m_refPos[x]==refPos) return x;
    return -1;
}


/*

 DINDELSEQUENCEHASH
 
 */

DindelSequenceHash::DindelSequenceHash(const std::string & sequence)
{
    makeHash(sequence);
}

void DindelSequenceHash::print() const
{
    for (Hash::const_iterator it=m_hash.begin();it!=m_hash.end();it++)
    {
        std::cerr << " hash: " << it->first << " => ";
                for (std::list<int>::const_iterator i=it->second.begin();i!=it->second.end();i++) std::cerr << " " << *i; std::cerr << std::endl;
    }
    
}

/*
int DindelSequenceHash::lookup(unsigned int key) const
{
    UniqueHash::const_iterator iter = m_uniqueHash.find(key);
    if (iter!=m_uniqueHash.end()) return iter->second; else return -1;
}
*/

const std::list<int> & DindelSequenceHash::lookupList(unsigned int key) const
{
    Hash::const_iterator iter = m_hash.find(key);
    if (iter == m_hash.end()) return emptyList; else return iter->second;
}

void DindelSequenceHash::makeHash(const std::string & sequence)
{
    for (size_t x=0;x<sequence.size()-DINDEL_HASH_SIZE;x++) m_hash[convert(sequence,x)].push_back(x);
}


/*
 *
 * DINDELWINDOW
 *
 *
 */


DindelWindow::DindelWindow(const std::vector<DindelVariant> & variants, const std::string & refHap, int refHapStart)
{

    
  
    m_pHaplotype_ma = NULL;

    m_candHapAlgorithm = LINEAR;

    // get refseq and set reference haplotype
    setupRefSeq(refHap, refHapStart);


    // setup window
    makeWindow(variants);



    // generate candidate haplotypes
    if (DINDEL_DEBUG) std::cout << "Generate candidate haplotypes." << std::endl;
    generateCandidateHaplotypes(variants);

    if (DINDEL_DEBUG) std::cout << "DindelWindow::DindelWindow DONE" << std::endl;
}


DindelWindow::DindelWindow(const std::vector<std::string> & haplotypeSequences, const std::string & refHap, int refHapStart, const std::string & refName)
{


    // setup reference haplotype. Will be the first in the vector
    setupRefSeq(refHap, refHapStart);
    m_chrom = refName;
    m_pHaplotype_ma = NULL;

    // globally align haplotypes to the reference sequence
    std::vector< MAlignData > maVector;

    
    
    std::set<std::string> uniqueHaplotypes;
    for (size_t h=0;h<haplotypeSequences.size();h++) uniqueHaplotypes.insert(haplotypeSequences[h]);

    // JTS: This assertion trips, change to warning
    //assert(uniqueHaplotypes.size() == haplotypeSequences.size());
    if(uniqueHaplotypes.size() != haplotypeSequences.size())
        std::cerr << "Warning: non-unique haplotype found\n";

    for (size_t h=0;h<haplotypeSequences.size();h++)
    {
        MAlignData _ma;
        _ma.position = 0;
        _ma.str = haplotypeSequences[h];

        std::stringstream ss; ss << "haplotype-" << h+1;

        _ma.name = ss.str();
        _ma.expandedCigar = StdAlnTools::expandCigar(StdAlnTools::globalAlignmentCigar(haplotypeSequences[h], refHap));

        if (DINDEL_DEBUG) std::cout << "DindelWindow::DindelWindow globalAlignmentCigar " << h << " vs root: " << _ma.expandedCigar << std::endl;
        
        maVector.push_back(_ma);
    }

    m_pHaplotype_ma = new MultiAlignment(refHap, maVector, refName);

    MultiAlignment & ma = *m_pHaplotype_ma;

    // get row indices for candidate haplotypes.
    size_t numHaplotypes = haplotypeSequences.size();
    std::vector<size_t> rowIdx(numHaplotypes,0);
    for(size_t h=0;h<numHaplotypes;h++) rowIdx[h] = ma.getIdxByName(maVector[h].name);

    
    size_t refRow = ma.getRootIdx();
    size_t numCols = ma.getNumColumns();
    
    assert(numCols>0);



     // Iterate over every column of the multiple alignment and detect variants

    for(size_t h = 0;h != numHaplotypes; h++)
    {
        size_t varRow = rowIdx[h];
        if (DINDEL_DEBUG) std::cout << "DindelWindow::DindelWindow Adding haplotype " << h << " in varRow " << varRow << std::endl;
        addHaplotype(DindelHaplotype(m_chrom, m_refSeq, refHapStart, *m_pHaplotype_ma, varRow, refRow));
    }
        
    delete m_pHaplotype_ma;
    m_pHaplotype_ma = NULL;


}

DindelWindow::~DindelWindow()
{
    if (m_pHaplotype_ma!=NULL) delete m_pHaplotype_ma;
}

void DindelWindow::generateCandidateHaplotypesFromMultiAlignment()
{
    

    


}


void DindelWindow::makeWindow(const std::vector<DindelVariant>& variants)
{

    if (variants.empty()) throw std::string("DindelWindow::no_variants");

    std::string chrom = "-1";
    int maxPos=4000000000, minPos=-1;

    for (size_t v=0;v<variants.size();v++)
    {
        DindelVariant var = variants[v];
        std::string vchrom = var.getChrom();

        // delOffs makes sure that the part being deleted by the variant is included in the window.
        int delOffs = 0;
        if (var.getLengthDifference()<0) delOffs = -var.getLengthDifference();

        int min_vpos = var.getPos();
        int max_vpos = min_vpos+delOffs;


        if (v==0)
        {
            maxPos = max_vpos;
            minPos = min_vpos;
            chrom = vchrom;
        } else
        {
            if (chrom.compare(vchrom)!=0) throw std::string("DindelWindow::multiple_chromosomes_in_same_window");
            if (max_vpos>maxPos) maxPos=max_vpos;
            if (min_vpos<minPos) minPos=min_vpos;
        }

        //TODO Need to check if variant is unique and it would be nice to add realignment of candidates here.
     }

    assert(minPos>m_leftRef);
    assert(maxPos<m_rightRef);

    m_chrom = chrom;

    if (DINDEL_DEBUG) 
    {
        std::cout << "window: chrom " << m_chrom << " left: " << m_leftRef << " right: " << m_rightRef << std::endl;
    }        
    if (1) std::cerr << "window size: " << m_rightRef-m_leftRef << std::endl;

}

void DindelWindow::setupRefSeq(const std::string & refHap, int refHapStart)
{
    
    m_leftRef = refHapStart;
    m_rightRef = refHapStart+int(refHap.size())-1;
    assert(m_leftRef<=m_rightRef);
    m_refSeq = refHap;

    
    haplotypes.clear();

    // first add reference haplotype
    haplotypes.push_back(DindelHaplotype(m_refSeq, m_leftRef, true));


}

DindelHaplotype DindelWindow::makeAltHapFromReference(const DindelVariant & var)
{
    int relPos = var.getPos() - m_leftRef;
    if (relPos<0)
    {
        std::cerr << "relPos: " << relPos << " var.getPos(): " << var.getPos() << " m_leftRef: " << m_leftRef << std::endl;
        assert("DindelWindow::generate_candidate_haplotypes_window_error");
    }
    // check if reference part of variant is consistent with its position
    if (m_refSeq.compare(relPos, var.getRef().length(), var.getRef())!=0)
    {
        throw std::string("DindelWindow::variant_refandpos_inconsistent_with_reference");
    }
    DindelHaplotype altHap(m_refSeq, m_leftRef, true);

    // note that addVariant may decide to not add the variant, for instance, when it conflicts with previously added variants
    altHap.addVariant(var);

    return altHap;
}

bool DindelWindow::addHaplotype(const DindelHaplotype & haplotype)
{
    if (m_hashAltHaps.find(haplotype.getSequence())==m_hashAltHaps.end())
    {
        haplotypes.push_back(haplotype);
        if (DINDEL_DEBUG)
        {
            std::cerr << "DindelWindow::addHaplotype candidate_haplotype  " << haplotypes.size()-1 << std::endl;
            haplotype.write(std::cerr);
            std::cerr << "DindelWindow::addHaplotype DONE" << std::endl;
        }
        m_hashAltHaps[haplotype.getSequence()] = haplotypes.size()-1;

        return true;
    }
    return false;
}

void DindelWindow::addVariant(const DindelVariant& variant, bool addToAll)
{
    if (DINDEL_DEBUG) std::cerr << "DindelWindow::addVariant START" << std::endl;
    
    // makeAltHapFromReference checks the variant against the reference sequence
    DindelHaplotype toRef = makeAltHapFromReference(variant);
    DindelVariant var=variant;
    bool added=addHaplotype(toRef);
    if (added) 
    {
            if (DINDEL_DEBUG) std::cout << "Adding variant " << var.getID() << std::endl;

            
            if (addToAll)
            {   // add to haplotypes other than the reference haplotype;
                int numHaps = int(haplotypes.size())-1;
                for (int h=1;h<numHaps;h++)
                {
                    // this DindelHaplotype constructor copies the m_refPos and m_seq etc from the haplotype
                    DindelHaplotype altHap(haplotypes[h],0);
                    if (altHap.addVariant(var)) 
                    {
                        addHaplotype(altHap); // only add haplotype if the addVariant was successful                        
                    }
                }
            }
    }
    if (DINDEL_DEBUG) std::cerr << "DindelWindow::addVariant END" << std::endl;

}

void DindelWindow::generateCandidateHaplotypes(const std::vector<DindelVariant> & variants)
{
    if (DINDEL_DEBUG) 
    {
        std::cout << "reference_haplotype: " << std::endl;
        haplotypes[0].write(std::cout);
    }


    if (m_candHapAlgorithm == LINEAR)
    {
        if (DINDEL_DEBUG) std::cout << "m_candHapAlgorithm LINEAR" << std::endl;
        // make one candidate haplotype per candidate variant
        for (size_t v=0;v<variants.size();v++)
        {
            addVariant(variants[v], false);
        }
    
    } else
    {
        throw std::string("DindelWindow::generateCandidateHaplotypes_unknown_algorithm");
    }
}


void DindelWindow::filterHaplotypes(const std::vector<bool> & filterHaplotype)
{
    assert (filterHaplotype.size() == haplotypes.size());
    std::vector<DindelHaplotype> newHaplotypes;
    for (size_t h=0;h<filterHaplotype.size();h++)
    {
        if (filterHaplotype[h])
        {
            // erase it from the hash

            std::map<std::string, int>::iterator it = m_hashAltHaps.find(haplotypes[h].getSequence());
            assert(it!=m_hashAltHaps.end());
            m_hashAltHaps.erase(it);

        } else {
            newHaplotypes.push_back(haplotypes[h]);
        }
    }
    haplotypes.swap(newHaplotypes);
}


/*


    DindelRealignWindowResult



 */

// Need to add pointer to VCF Header instance
void DindelRealignWindowResult::Inference::outputAsVCF(const DindelVariant & var, const DindelRealignWindowResult & result, std::ostream& out) const
{
    int iqual = (qual<0.0)?0:int(qual);
    std::string filter="NoCall";
    if (iqual>10 && qual<20)
        filter = "LowQuality";
    else if (iqual>=20)
        filter = "PASS";

    std::set<int> hps;
    hps.insert(result.m_pDindelRealignWindow->getDindelWindow().getHaplotypes()[0].getHomopolymerLengthRefPos(var.getPos()+1));
    hps.insert(result.m_pDindelRealignWindow->getDindelWindow().getHaplotypes()[0].getHomopolymerLengthRefPos(var.getPos()+0));
    hps.insert(result.m_pDindelRealignWindow->getDindelWindow().getHaplotypes()[0].getHomopolymerLengthRefPos(var.getPos()-1));
    int hp = *hps.rbegin();

    out.precision(5);
    out.setf(std::ios::fixed,std::ios::floatfield);
    out << var.getChrom() << "\t" << var.getPos() << "\t.\t" << var.getRef() << "\t" << var.getAlt() << "\t" << ( iqual ) << "\t" << filter << "\t";
    out << "AF=" << freq;
    out << ";DP=" << numReadsForward+numReadsReverse << ";NF=" << numReadsForward << ";NR=" << numReadsReverse << ";NF0=" << this->numReadsForwardZeroMismatch;
    out << ";NR0=" << this->numReadsReverseZeroMismatch << ";NU=" << numUnmapped << ";NH=" << result.numHaplotypes << ";HSR=" << result.numHapSpecificReads;

    out << ";HPLen=" << hp;
    out << ";SB=" << this->strandBias;
    out << ";NumLib=" << this->numLibraries;
    out << ";NumFragments=" << this->numReadNames;
    out << ";SingleSample=" << (this->isSingleSampleCall?"1":"0");
    out << ";SingleSampleQual=" << int(this->singleSampleQual);

    if (!infoStr.empty()) out << ";" << this->infoStr;

    if (!this->histDistance.empty())
    {
        out << ";HistDist=";
        for (size_t x=0;x<histDistance.size();x++)
        {
            if (x>0) out << ",";
            out << histDistance[x];
        }
    }
    
    if (!this->histAlignLik.empty())
    {
        out << ";HistLik=";
        for (size_t x=0;x<histAlignLik.size();x++)
        {
            if (x>0) out << ",";
            out << histAlignLik[x];
        }
    }

    if (!this->histMapQ.empty())
    {
        out << ";HistMapQ=";
        for (size_t x=0;x<histMapQ.size();x++)
        {
            if (x>0) out << ",";
            out << histMapQ[x];
        }
    }

   // output genotyping information

    if (result.outputGenotypes)
    {
        // JS: This has been hacked out since we now output to an ostream instead of a VCFFile
        assert(false);
#if 0
        out << "\tGT:GQ:GL";
        const std::vector<std::string> & samples = vcfFile.getSamples();
        for (size_t x=0;x<samples.size();x++)
        {
            const std::string & sample = samples[x];
        const  DindelRealignWindowResult::SampleToGenotypes::const_iterator sample_it=result.sampleToGenotypes.find(sample);
        out << "\t";
            if (sample_it!=result.sampleToGenotypes.end())
            {
        DindelRealignWindowResult::VarToGenotypeCall::const_iterator it_var = sample_it->second.find(var);
                assert(it_var!=sample_it->second.end());
                const DindelRealignWindowResult::GenotypeCall & gc = it_var->second;

                if (gc.count == 0) out << "0/0";
                else if (gc.count == 1) out << "0/1";
                else if (gc.count == 2) out << "1/1";
                out << ":" << int(gc.qual);
        //out <<  std::iostream::fixed;
                out << std::setprecision(5) << ":" << gc.gl[0] << "," << gc.gl[1] << "," << gc.gl[2];
            } else
            {
                out << "."; // no call
            }
        }
#endif
    }
    out << std::endl;
}

double DindelRealignWindowResult::Inference::computeStrandBias(int numForward, int numReverse)
{
    // computes Bayes factor for strand bias
    // Even in the case of no strand bias, the expected frequency is allowed to differ from 0.5 in order to allow for
    // mapper bias.
    int total = numForward+numReverse;
    
    const int N=50;
    double p = (1.0)/double(N);

    double bBias = 0.0;
    double bNoBias = 0.0;

    double nf = numForward;
    double nr = numReverse;
    double comb=(total<=10) ? double(combinations(total, numForward )) : exp( lgamma(double(total+1))-lgamma(nf+1.0)-lgamma(double(total-numForward+1)));

    double tBias = 0.35; // point at which we believe there is a bias

    for (int i=0;i<N/2;i++)
    {
        double q = double(i+1)*tBias/(double(N+1)/2.0);
        bBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;

        q = 1.0-double(i+1)*tBias/(double(N+1)*2.0);
        bBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
    }

    
    for (int i=0;i<N/2;i++)
    {
        double q = 0.5+(0.5-tBias)*double(i+1)/(double(N+1)/2.0);
        bNoBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;

        q = 0.5-(0.5-tBias)*double(i+1)/(double(N+1)/2.0);
        bNoBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
    }
    
    
    double bf = log(bBias)-log(bNoBias);
    return bf;
}

void DindelRealignWindowResult::Inference::addDistanceToHistogram(int distance)
{
    const int D=25;
    if (this->histDistance.empty()) histDistance = std::vector<int>(D,0);
    if (distance>=0) {
        if (distance>D-1) distance = D-1;
        histDistance[distance]++;
    }
}

/// logLik is assumed to be the average per base
void DindelRealignWindowResult::Inference::addAlignLikToHistogram(double logLik)
{
    const int D=25;
    if (this->histAlignLik.empty()) histAlignLik = std::vector<int>(D,0);
    int bin = int(round(logLik* 25.0/(-10.0/50.0)));
    assert(bin>=0);
    if (bin>D-1) bin=D-1;
    histAlignLik[bin]++;
}

void DindelRealignWindowResult::Inference::addMapQToHistogram(double mappingQuality)
{
    const int D=8;
    if (this->histMapQ.empty()) histMapQ = std::vector<int>(D,0);
    int bin = int(round(mappingQuality/10.0)); 
    if (bin>D-1) bin=D-1;
    if (bin>=0) histMapQ[bin]++;
}

void DindelRealignWindowResult::outputVCF(std::ostream& out)
{
    VarToInference::const_iterator iter = variantInference.begin();

    for (;iter!=variantInference.end();iter++)
    {
        iter->second.outputAsVCF(iter->first, *this, out);
    }    
}


DindelRealignWindow::ReadSamples::ReadSamples(const std::vector<DindelRead>& dindelReads)
{

    typedef HashMap<std::string, std::list<int> > Sample_to_read;
    Sample_to_read sample_to_read;
    for (int r=0;r<int(dindelReads.size());r++)
    {
        const std::string & sample = dindelReads[r].getSampleName();
        sample_to_read[sample].push_back(r);
    }

    readIndices.reserve(dindelReads.size()+sample_to_read.size());
    int sample_idx = 0;
    for (Sample_to_read::const_iterator it = sample_to_read.begin();it!=sample_to_read.end();it++, sample_idx++)
    {
        samples.push_back(it->first);        
        for (std::list<int>::const_iterator rit = it->second.begin();rit!=it->second.end();rit++)
        {
            readIndices.push_back(*rit);
        }
        readIndices.push_back(-1);
    }

    int numSamples = std::max(1, int(samples.size()));
    // compute number of states
    std::vector<double> Q1(numSamples*2+1,0.0),Q2(numSamples*2+1,0.0); 
    const double lc = log(2.0);

    Q1[0] = 0.0;
    Q1[1] = 0.0+lc;
    Q1[2] = 0.0;
    
    for (int s=2;s<=numSamples;s++)
    {

       Q2[0] = Q1[0]+0;
       Q2[1] = addLogs(Q1[1]+0,lc+Q1[0]+0);
       for (int k=2;k<=2*s-2;k++)
       {
        Q2[k] = addLogs(addLogs(Q1[k]+0,lc+Q1[k-1]+0),Q1[k-2]+0);
       }
       Q2[2*s-1] = addLogs(lc+Q1[2*s-1-1]+0,Q1[2*s-1-2]+0);
       Q2[2*s] = Q1[2*s-2]+0;
       Q1.swap(Q2);
    }
    logNumStates = Q1;

}



/*
 *
 * DINDELREALIGNWINDOW
 *
 *
 */

DindelRealignWindow::DindelRealignWindow(const DindelWindow* pDindelWindow,
        std::vector<DindelRead> & dindelReads,
        const DindelRealignParameters & dindelRealignParameters) :
    m_dindelWindow(*pDindelWindow),
    m_pDindelReads(&dindelReads),
    realignParameters(dindelRealignParameters),
    m_readSamples(dindelReads)
{ 
    if (DINDEL_DEBUG) std::cout << "DindelRealignWindow::DindelRealignWindow STARTED" << std::endl;
    if (DINDEL_DEBUG) std::cout << "DindelRealignWindow::DindelRealignWindow DONE" << std::endl;    
}

void DindelRealignWindow::run(const std::string & algorithm, std::ostream& out)
{
    if (algorithm == "hmm")
    {
        algorithm_hmm(out);
    }
    else
    {
        throw std::string("DindelRealignWindow::run::unknown_algorithm");
    }

}

void DindelRealignWindow::addSNPsToHaplotypes()
{
    typedef std::map<size_t, std::vector<DindelVariant> > FreqCandidates;
    FreqCandidates freqCandidates;

    for (PosToCandidateSNP::const_iterator iter = m_posToCandidateSNP.begin(); iter != m_posToCandidateSNP.end(); iter++)
    {
        for (size_t baseIdx=0;baseIdx<4;baseIdx++)
        {
            if (iter->second.baseToReads[baseIdx].size()>=realignParameters.addSNPMinCount)
            {
                std::string ref; ref+=iter->second.m_refBase;
        std::string alt; alt+=iter->second.getBase(baseIdx);
        DindelVariant snp(m_dindelWindow.getChrom(), ref, alt, iter->first);
        snp.setPriorProb(realignParameters.variantPriors.getDefaultProbVariant(snp.getType()));
                freqCandidates[iter->second.baseToReads[baseIdx].size()].push_back(snp);
        if (!QUIET) std::cerr << "DindelRealignWindow::addSNPsToHaplotypes adding SNP " << m_dindelWindow.getChrom() << " " << ref << " "  <<  alt << " " <<  iter->first << std::endl;
            }
        }
    }

    // add snps until maximum is reached.
    int numAdded=0;
    for (FreqCandidates::const_reverse_iterator iter=freqCandidates.rbegin();iter!=freqCandidates.rend();iter++)
    {
        bool stop=false;
        for(size_t x=0;x<iter->second.size();x++)
        {
            m_dindelWindow.addVariant(iter->second[x], true);
            numAdded++;
            if (numAdded == realignParameters.addSNPMaxSNPs)
            {
                stop =true;
                break;
            }

        }
        if (stop) break;
    }

}



void DindelRealignWindow::algorithm_hmm(std::ostream& out)
{
    if (DINDEL_DEBUG) std::cerr << "DindelRealignWindow::algorithm_hmm STARTED" << std::endl;

    // filter the haplotypes
    filterHaplotypes();


    size_t numHaps = m_dindelWindow.getHaplotypes().size();
    assert(numHaps>=1);

    // process the user-defined haplotypes
    computeReadHaplotypeAlignmentsUsingHMM(0, numHaps-1);

    // now add SNPs identified by the high quality ungapped alignments
    if (ADDSNPS && realignParameters.addSNPMaxSNPs>0)
    {
        addSNPsToHaplotypes();

        // compute the likelihoods for the additional haplotypes
        computeReadHaplotypeAlignmentsUsingHMM(numHaps, m_dindelWindow.getHaplotypes().size()-1);
    }

    // estimate haplotype frequencies using user specified max mapping thresholds
    DindelRealignWindowResult result;
    if (realignParameters.doEM)
        result = estimateHaplotypeFrequenciesModelSelection(realignParameters.minLogLikAlignToRef, realignParameters.minLogLikAlignToAlt, true, true);
    else
        result = estimateHaplotypeFrequencies(realignParameters.minLogLikAlignToRef, realignParameters.minLogLikAlignToAlt, true, true);
    // estimate variant qualities using the standard q40 and q60 thresholds
    int maxMapQuals[] = { 25, 80 };
    if (0) for (int m=0;m<2;m++)
    {
        double minLogLikAlign=-((double) maxMapQuals[m])*.23026-log(double(DINDEL_HMM_BANDWIDTH));
            bool print = false;
            if (m==0) print = true;
        DindelRealignWindowResult result_q40 = estimateHaplotypeFrequencies(minLogLikAlign, minLogLikAlign, false, print);
        // update info strings in result
        for (DindelRealignWindowResult::VarToInference::const_iterator it = result_q40.variantInference.begin(); it!= result_q40.variantInference.end(); it++)
        {
        DindelRealignWindowResult::VarToInference::iterator it_def = result.variantInference.find(it->first);
        assert(it_def!=result.variantInference.end());
        std::ostringstream os; os << "QMQ" << maxMapQuals[m] << "=" << (int(it->second.qual)<0?0:int(it->second.qual));
        os << ";DPMQ" << maxMapQuals[m] << "=" << int(it->second.numReadsForward+it->second.numReadsReverse);
        if (it_def->second.infoStr.empty()) 
            it_def->second.infoStr = os.str();
        else
            it_def->second.infoStr += ";" + os.str();

        }
    }

    if (realignParameters.genotyping)
    {
        addDiploidGenotypes(result, true);
        result.outputGenotypes=true;
    } else
    {
        result.outputGenotypes=false;
    }


    // output the results. result does marginalization over haplotypes
    result.outputVCF(out);

    if (DINDEL_DEBUG) std::cerr << "DindelRealignWindow::algorithm_hmm DONE" << std::endl;

}



void PosToCandidateSNP::addSNP(int readIndex, int refPos, char refBase, char altBase)
{

   //std::cerr << " PosToCandidateSNP::size(): " << this->size() << " readIndex: " << readIndex << " refPos: " << refPos << " refBase: " << refBase << " altBase: " << altBase << std::endl;

    PosToCandidateSNP::iterator iter = this->find(refPos);
    if (iter == this->end())
    {
        std::pair<PosToCandidateSNP::iterator, bool> pib  = this->insert(PosToCandidateSNP::value_type(refPos, CandidateSNP(refBase)));
        iter = pib.first;
    }

    iter->second.addRead(refBase, altBase, readIndex);

}

void DindelRealignWindow::filterHaplotypes()
{
   
}



ReadHaplotypeAlignment DindelRealignWindow::computeReadHaplotypeAlignment(size_t readIndex, const DindelHaplotype& haplotype, int start, int end, const std::vector<double> & lpError, const std::vector<double> & lpCorrect, bool rcReadSeq)
{
    const DindelRead & read = (*m_pDindelReads)[readIndex];
    // compute the log-likelihood of a single candidate alignment

    const int DRHA=0;

    if (DINDEL_DEBUG &&DRHA) std::cerr << std::endl <<  "====> START 1 DindelRealignWindow::computeReadHaplotypeAlignment " << read.getID() << std::endl;
    // scenarios
    /* Capital R's are matches, r's are mismatches
    1. The seeds are too long to detect the possible gapped alignment of read to haplotype
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                     rrrrrrrRRRRRRRRRRRRRRRRRRRRRRRRRRR
        which should be
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                  RRR   rrrrRRRRRRRRRRRRRRRRRRRRRRRRRRR


    2. Deletion in read
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                   RRRRRRR      RRRRRRRRRRRRRRRRRRRRRRRR
    */

    // start and end are inclusive
    // start may be <0, and end may be > haplotype length
    // lpError, lpCorrect

    int dlen = end+1-start-read.length();
    int rlen = read.length();
    int hlen = haplotype.length();
    const double minLpCorrect = realignParameters.addSNPMinLogBaseQual;
   
    if (DINDEL_DEBUG && DRHA)
    {
        std::cerr << std::endl;
        std::cerr << "HAPLOTYPE: " << haplotype.getSequence() << " isREF " << haplotype.isReference() << std::endl;
        std::string pad;
        std::string rseq = read.getSequence();
        if (start>=0) 
        {
            pad = std::string(start, ' ');
        }
        else {
           rseq = read.getSequence().substr(-start, read.length());
        }       

        std::cerr << "READ:      " << pad << rseq << std::endl;
        std::cerr << "DLEN: " << dlen << " RLEN: " << rlen << " hlen: " << hlen << " start: " << start << " end: " << end << " read_index: " << readIndex << " " << read.getID() <<  std::endl;
    }

    int nmm=-1; // number of mismatches between read and candidate haplotype
    double logLik;
    bool isUngapped = false;
    if (dlen == 0)
    {
        // possible ungapped alignment
        typedef std::pair<int, char> Mismatch;
        std::list<Mismatch> mismatches;
        bool countMismatch = (realignParameters.addSNPMinMappingQual<read.getMappingQual())?true:false;

        double ll=0.0;
        int b=0;
    double all_match = 0.0;
    nmm=0;
        if (!rcReadSeq) { 
        for (int x=start;x<=end;x++,b++)
            {
        bool match = true;
        if (x>=0 && x<=hlen-1) 
        {
            match = (read.getBase(b) == haplotype.getSequence()[x])?true:false;
        }
        if (match) ll += lpCorrect[b]; else { 
            if (DINDEL_DEBUG) std::cerr << " mm " << b << " lpError: " << lpError[b] << std::endl;
                    ll += lpError[b]; nmm++;
                    if (countMismatch && lpCorrect[b]>minLpCorrect) mismatches.push_back(Mismatch(x, read.getBase(b)));
                }
        all_match += lpCorrect[b];
        }
    } else {
        for (int x=start;x<=end;x++,b++)
        {
            bool match = true;
            if (x>=0 && x<=hlen-1)
        {
            match = (_complement(read.getBase(rlen-1-b)) == haplotype.getSequence()[x])?true:false;
        }
            if (match) ll += lpCorrect[rlen-1-b]; else { 
                    ll += lpError[rlen-1-b]; nmm++;
                    if (countMismatch && lpCorrect[rlen-b-1]>minLpCorrect) mismatches.push_back(Mismatch(x, _complement(read.getBase(rlen-1-b))));
                }
        }
    }

        if (nmm<=realignParameters.addSNPMaxMismatches) {
            for (std::list<Mismatch>::const_iterator iter = mismatches.begin(); iter != mismatches.end(); iter++) {
                int refPos = haplotype.getRefBase(iter->first);
        if (refPos>0) 
        {
                    m_posToCandidateSNP.addSNP(readIndex, refPos, haplotype.getSequence()[iter->first], iter->second); // pos, ref, alt resp.
        }
            }
        }
    logLik = ll;
    isUngapped = true;
    }
    else { assert(1==0); }


    assert(!isnan(logLik));
    assert(logLik<0.0);

    logLik -= log(double(DINDEL_HMM_BANDWIDTH)); // added for consistency with the HMM, which has a prior 

    // cap likelihood based on mapping quality
    // addLogs(logLik, -read.getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));

    if (DINDEL_DEBUG && DRHA)
    {
    std::cerr << "LOGLIK: " << logLik << std::endl;
    }
    if (DINDEL_DEBUG &&DRHA) std::cerr << std::endl <<  "====> END DindelRealignWindow::computeReadHaplotypeAlignment" << std::endl;
    return ReadHaplotypeAlignment(logLik, nmm, isUngapped);
}



void DindelRealignWindow::HMMAlignReadAgainstHaplotypes(size_t readIndex, size_t firstHap, size_t lastHap, const std::vector<double> & lpCorrect, const std::vector<double> & lpError)
{

    const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();

    DindelRead & read = *(m_pDindelReads->begin()+readIndex);

    for (size_t h=firstHap;h<=lastHap;++h)
    {
        const DindelHaplotype & haplotype = haplotypes[h];
        DindelHMM hmm(read, haplotype);
        ReadHaplotypeAlignment rha_hmm = hmm.getAlignment();
    if (DINDEL_DEBUG) std::cerr << " HMM read " << readIndex << " haplotype " <<  h  << " loglik: " << rha_hmm.logLik << " postProb: " << rha_hmm.postProbLastReadBase << " hapPosLastReadBase: " << rha_hmm.hapPosLastReadBase << std::endl;

        // attempt ungapped alignment for detecting SNPs
        if (rha_hmm.postProbLastReadBase>this->realignParameters.minPostProbLastReadBaseForUngapped)
        {
    

            int end = rha_hmm.hapPosLastReadBase;
            int start = end-read.length()+1;
            bool rcReadSeq=read.getRCRead();
            ReadHaplotypeAlignment rha_ung = computeReadHaplotypeAlignment(readIndex, haplotype, start, end, lpError, lpCorrect, rcReadSeq);

            if (rha_ung.nmm<=2)
            {
                rha_hmm.isUngapped=true;
                rha_hmm.nmm = rha_ung.nmm;
            }
        if (DINDEL_DEBUG) std::cerr << " HMM read " << readIndex << " haplotype " << h  << " ungapped loglik: " << rha_ung.logLik << " nmm: " << rha_ung.nmm << std::endl;
        }
    
        hapReadAlignments[h][readIndex]=rha_hmm;
    }

}

void DindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM(size_t firstHap, size_t lastHap)
{

    if (DINDEL_DEBUG) std::cerr << "DindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM firstHap " << firstHap << " " << lastHap << std::endl;
    // store information whether reads were realigned against all haplotypes or just the reference haplotype?
    const std::vector<DindelHaplotype> haplotypes = m_dindelWindow.getHaplotypes();

    if (m_readMapsToReference.empty()) m_readMapsToReference = std::vector<int>(m_pDindelReads->size(),0);

    assert(haplotypes.size()>lastHap);
    for (size_t h=firstHap;h<=lastHap;h++)
    {
        if (haplotypes[h].isReference())
            hapReadAlignments.push_back( std::vector<ReadHaplotypeAlignment>(m_pDindelReads->size(), ReadHaplotypeAlignment(realignParameters.minLogLikAlignToRef,-1)));
        else
        hapReadAlignments.push_back( std::vector<ReadHaplotypeAlignment>(m_pDindelReads->size(), ReadHaplotypeAlignment(realignParameters.minLogLikAlignToAlt,-1)));
    }

    for (size_t r=0;r<m_pDindelReads->size();r++) if (m_readMapsToReference[r]==0)
    {
    const DindelRead & read = (*m_pDindelReads)[r];
    if  (DINDEL_DEBUG) std::cerr << "\n*****\nDindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM reads[" << r << "]: " << read.getID() << std::endl;
    
        // only realign if it has not yet been found to map to the reference
        std::vector<double> lpCorrect, lpError;
        read.getLogProbCorrectError(lpCorrect, lpError);

        HMMAlignReadAgainstHaplotypes(r, firstHap, lastHap, lpCorrect, lpError);
    }
}

double DindelRealignWindow::getHaplotypePrior(const DindelHaplotype & h1, const DindelHaplotype & h2) const
{
   const std::vector<DindelVariant> & v1 = h1.getVariants();
   const std::vector<DindelVariant> & v2 = h2.getVariants();
   std::set<DindelVariant> uniqueVariants;
   for (size_t x=0;x<v1.size();x++) uniqueVariants.insert(v1[x]);
   for (size_t x=0;x<v2.size();x++) uniqueVariants.insert(v2[x]);

   double logPrior = 0.0;
  
   for (std::set<DindelVariant>::const_iterator it1=uniqueVariants.begin();it1!=uniqueVariants.end();it1++)
       logPrior += log(it1->getPriorProb());

   return logPrior;  
}

void DindelRealignWindow::addDiploidGenotypes(DindelRealignWindowResult& result, bool useEstimatedHaplotypeFrequencies = false)
{
   const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   const std::vector<DindelRead> & reads = *m_pDindelReads;
   if (DINDEL_DEBUG) std::cerr << "GENOTYPING START" << std::endl;
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::addDiploidGenotypes incomplete haplotype alignments.");
   int numHaps = int(haplotypes.size());
   
   typedef HashMap<std::string, std::list<int> > SampleToReads;
   SampleToReads sampleToReads;

   for (size_t r=0;r<reads.size();r++)
   {
       const std::string & sampleName = reads[r].getSampleName();
       sampleToReads[sampleName].push_back(r);
   }

   std::vector<int> allowedHaplotype(haplotypes.size(),1);
   useEstimatedHaplotypeFrequencies = false;
   /*
   if (useEstimatedHaplotypeFrequencies)
   {
       // remove haplotypes with frequency below threshold
       if (result.haplotypeFrequencies.size()!=haplotypes.size()) 
       {
           std::cerr << "result.haplotypeFrequencies.size: " << result.haplotypeFrequencies.size() << " haplotypes.size: " << haplotypes.size() << std::endl;        
           throw std::string("frequencies and haplotypes inconsistent.");
       }
       for (size_t h=1;h<haplotypes.size();h++) if (result.haplotypeFrequencies[h]<1e-8) allowedHaplotype[h]=0;
   }
   */

   typedef std::map<DindelVariant, HashMap<int, int> > VariantToHaplotype;
   VariantToHaplotype variantToHaplotype;

   for (int i=0;i<numHaps;i++)
   {
       const std::vector<DindelVariant> & vars = haplotypes[i].getVariants();
       for (size_t x=0;x<vars.size();x++) variantToHaplotype[vars[x]][i]=1;
   }
   

   for (SampleToReads::const_iterator sample_it = sampleToReads.begin(); sample_it != sampleToReads.end(); sample_it++)
   {
       const std::string & sample = sample_it->first;

       // infer likelihoods for each pair of haplotypes
       std::vector<double> llHapPairs(numHaps*numHaps);

       double max_ll=-HUGE_VAL;
       int h1max=-1, h2max=-1;

       for (int h1=0;h1<numHaps;h1++) if (allowedHaplotype[h1])
           for (int h2=h1;h2<numHaps;h2++) if (allowedHaplotype[h2])
           {
               double ll=0.0;
               for (std::list<int>::const_iterator _r = sampleToReads[sample].begin(); _r != sampleToReads[sample].end(); ++_r)
               {
                   int r = *_r;
                   double mq = -reads[r].getMappingQual()*0.23056-log(double(DINDEL_HMM_BANDWIDTH));
                   double lh1 = addLogs(hapReadAlignments[h1][r].logLik, mq);
                   double lh2 = addLogs(hapReadAlignments[h2][r].logLik, mq);
                   ll += log(.5) + addLogs(lh1, lh2);
               }

               /*
               if (!useEstimatedHaplotypeFrequencies)
               {
                    ll += getHaplotypePrior(haplotypes[h1], haplotypes[h2]);
               } else 
               {
                   double f1 = result.haplotypeFrequencies[h1]; if (f1<1e-8) f1=1e-8;
                   double f2 = result.haplotypeFrequencies[h2]; if (f2<1e-8) f2=1e-8;               
                   ll += (log(f1)+log(f2));
               }
               */


               llHapPairs[h1*numHaps+h2]=ll;
               llHapPairs[h2*numHaps+h1]=ll;

               if (ll>max_ll)
               {
                   max_ll = ll;
                   h1max = h1;
                   h2max = h2;
               }
           }

       assert(h1max!=-1 && h2max !=-1);
       // set genotypes based on best haplotype pair
       const std::vector<DindelVariant> & v1 = haplotypes[h1max].getVariants();
       const std::vector<DindelVariant> & v2 = haplotypes[h2max].getVariants();

       // genotype call for most likely pair of haplotypes
       result.sampleToGenotypes[sample] = DindelRealignWindowResult::VarToGenotypeCall();
       DindelRealignWindowResult::VarToGenotypeCall & gc = result.sampleToGenotypes[sample];

       // initialize to NOT called
       for (VariantToHaplotype::const_iterator it = variantToHaplotype.begin();it!=variantToHaplotype.end();it++)
       {
           gc[it->first] = DindelRealignWindowResult::GenotypeCall();
       }

       // gc contains genotypes for maximum likelihood pair of haplotypes
       for (size_t x=0;x<v1.size();x++) gc[v1[x]].count += 1;
       for (size_t x=0;x<v2.size();x++) gc[v2[x]].count += 1;

       DindelRealignWindowResult::VarToGenotypeCall::iterator it1;
       for (it1=gc.begin();it1!=gc.end();it1++)
       {
           // initialize quality
           it1->second.qual = 123456.0;
           it1->second.called = true;
           for (int x=0;x<3;x++) it1->second.gl[x] = -HUGE_VAL;
       }

       for (int h1=0;h1<numHaps;h1++) if (allowedHaplotype[h1])
           for (int h2=h1;h2<numHaps;h2++) if (allowedHaplotype[h2])
           {
               double ll_this_pair = llHapPairs[h1*numHaps + h2];

               if (true)// (h1!=h1max || h2!=h2max)) NOT NECESSARY?
               {
                   // different pair of haplotypes
                   // only consider likelihood difference if it is higher, otherwise genotype quality will be zero.
                    const std::vector<DindelVariant> & _v1 = haplotypes[h1].getVariants();
                    const std::vector<DindelVariant> & _v2 = haplotypes[h2].getVariants();
                    DindelRealignWindowResult::VarToGenotypeCall _gc;

                    for (size_t x=0;x<_v1.size();x++) _gc[_v1[x]].count += 1;
                    for (size_t x=0;x<_v2.size();x++) _gc[_v2[x]].count += 1;

                    // update genotype qualities for called variants
                    for (it1 = gc.begin();it1 != gc.end();it1++)
                    {
                        DindelRealignWindowResult::VarToGenotypeCall::const_iterator it2;
                        const DindelVariant & variant = it1->first;
                        it2 = _gc.find(variant);
                        if (it2 == _gc.end())
                        {
                            // variant not called in this pair of haplotypes
                            int count = 0;
                            it1->second.gl[count] = addLogs(it1->second.gl[count], ll_this_pair);
                        } else
                        {
                            it1->second.gl[it2->second.count] = addLogs(it1->second.gl[it2->second.count], ll_this_pair);
                        }

                    }
               }
           }


       for (it1 = gc.begin();it1 != gc.end();it1++)
       {
           int count = it1->second.count;
           double qual = 1000.0;
           for (int c=0;c<3;c++) 
           {
               double d = it1->second.gl[count]-it1->second.gl[c];
               if (c!=count && d<qual)
               {
                   qual = d;
               }
           }
           if (qual<0.0) qual = 0.0;
           it1->second.qual = 10.0*qual/log(10.0);
       }
   }
   if (DINDEL_DEBUG) std::cerr << "GENOTYPING END" << std::endl;

}

DindelRealignWindow::HapIdxToCoalescentResult DindelRealignWindow::computeHaplotypesCoalescent(std::list<int> & calledHaplotypes, std::vector<double> & calledFreqs, std::list<int> candidateHaplotypes, double minLogLik)
{
   const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   const std::vector<DindelRead> & reads = *m_pDindelReads;
   
   
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::addDiploidGenotypes incomplete haplotype alignments.");

   int numSamples = int (m_readSamples.samples.size());

   assert(calledHaplotypes.size()==calledFreqs.size());
   std::vector<double> logCalledFreqs(calledFreqs.size()+1,0.0); // the last one is for the candidate haplotype
   for (size_t x=0;x<calledFreqs.size();x++) 
   {
       logCalledFreqs[x]=log(calledFreqs[x]);
       //std::cerr << "calledFreqs[" << x << "]: " << calledFreqs[x] << std::endl;

   }

   /*
   typedef std::map<DindelVariant, HashMap<int, int> > VariantToHaplotype;
   VariantToHaplotype variantToHaplotype;

   for (int i=0;i<numHaps;i++)
   {
       const std::vector<DindelVariant> & vars = haplotypes[i].getVariants();
       for (size_t x=0;x<vars.size();x++) variantToHaplotype[vars[x]][i]=1;
   }
   */

   // setup prior
    
   double theta = this->realignParameters.variantPriors.m_probINDEL;

   std::vector<double> prior(2*numSamples+1, 0.0);

   double sum=0.0;
   for (int k=1;k<=2*numSamples;k++) sum += 1.0/double(k);
   sum *=theta;
   prior[0]=log((1.0-sum));
   //prior[2*numSamples]=log(.5*(1.0-sum));

   double varsum=-HUGE_VAL;
   for (int k=1;k<=2*numSamples;k++)
   {
       prior[k]=log(theta*(1.0/double(k)))-m_readSamples.logNumStates[k];
       varsum=addLogs(varsum,prior[k]);
   }
   double z = log(sum)-varsum;
   for (int k=1;k<=2*numSamples;k++) prior[k]-=z;

   double checkz = -HUGE_VAL;
   for (int k=0;k<=2*numSamples;k++) 
   {
       checkz = addLogs(checkz,prior[k]);
           if (DEBUG_CALLINDEL)
           {
                std::cerr << "prior[" << k << "]: " << prior[k] << std::endl;
           }
   }
   assert (fabs(checkz)<0.001);

   HapIdxToCoalescentResult result;

   std::vector<double> Q1(numSamples*2+1,0.0), Q2(numSamples*2+1,0.0);

   

   for (std::list<int>::const_iterator cand_it=candidateHaplotypes.begin();cand_it!=candidateHaplotypes.end();cand_it++)
   {
       
       double qual_sample_lik = -100.0;// quality of call using simple likelihood base threshold.


       // setup allowed haplotypes
       std::vector<int> allowedHaplotypes;
       for (std::list<int>::const_iterator it =calledHaplotypes.begin();it!=calledHaplotypes.end();it++) allowedHaplotypes.push_back(*it);
       allowedHaplotypes.push_back(*cand_it);

       int numAH = int(allowedHaplotypes.size());


       // compute genotype likelihoods

       std::vector<double> llHapPairs(numSamples*3,-HUGE_VAL);

       std::vector<double> tmpLL(numAH*numAH,0.0);

       int read_idx = 0;
       for (int sample_idx = 0; sample_idx < numSamples;sample_idx++)
       {
           // infer likelihoods for non-candidate haplotype/cand haplotype genotypes

           for (int x=0;x<numAH*numAH;x++) tmpLL[x]=0.0;
           while (true)
           {
               int r = m_readSamples.readIndices[read_idx];
               read_idx++;
               if (r==-1) break;
            // if (m_readMapsToReference[r]) continue; // we treat these reads as uninformative. if used bias to low frequencies, but will not result in undercalling
               double mq = addLogs(minLogLik,-reads[r].getMappingQual()*0.23056-log(double(DINDEL_HMM_BANDWIDTH)));

               for (int h1=0;h1<numAH;h1++)
                   for (int h2=h1;h2<numAH;h2++)
                   {                      
                       double lh1 = addLogs(hapReadAlignments[allowedHaplotypes[h1]][r].logLik, mq);
                       double lh2 = addLogs(hapReadAlignments[allowedHaplotypes[h2]][r].logLik, mq);
               //if (h1==0 && h2==1) {
            //       std::cerr << "r: " << r << " lh1: " << lh1 << " lh2: " << lh2 << std::endl;
               //}
                       double ll = log(.5) + addLogs(lh1, lh2);
                       tmpLL[h1*numAH+h2]+=ll;
                   }
               if (DEBUG_CALLINDEL) { 
               std::cerr << "AFTER READ " << r << std::endl;
               for (int x=0;x<numAH*numAH;x++) std::cerr << "    TMPLIK[" << x << "]:" << tmpLL[x] << std::endl;
           }
           }

           // integrate over the called haplotypes
           for (int h1=0;h1<numAH;h1++)
           {
               for (int h2=h1;h2<numAH;h2++)
               {
                   int gt = 0;
                   if (h1 == numAH-1) gt++;
                   if (h2 == numAH-1) gt++;

                   llHapPairs[sample_idx*3+gt] = addLogs(llHapPairs[sample_idx*3+gt], tmpLL[h1*numAH+h2]+logCalledFreqs[h1]+logCalledFreqs[h2]);
               }
           }
       
           // make call based on individual sample likelihoods
           if (llHapPairs[sample_idx*3+1]>llHapPairs[sample_idx*3+0]+this->realignParameters.singleSampleHetThreshold || llHapPairs[sample_idx*3+2]>llHapPairs[sample_idx*3+0]+this->realignParameters.singleSampleHomThreshold)
           {
               double d1=llHapPairs[sample_idx*3+1]-llHapPairs[sample_idx*3+0];
               double d2=llHapPairs[sample_idx*3+2]-llHapPairs[sample_idx*3+0];
               double q = (d1>d2)?d1:d2;
               if (qual_sample_lik<0.0) qual_sample_lik=0.0;
               qual_sample_lik += q;
           }

           if (DEBUG_CALLINDEL)
           {
                std::cerr << "haplotype[" << allowedHaplotypes[numAH-1] << "]: liklihoods " << sample_idx << " " << m_readSamples.samples[sample_idx] << " " << llHapPairs[sample_idx*3+0] << " " << llHapPairs[sample_idx*3+1] << " " << llHapPairs[sample_idx*3+2];
                if (llHapPairs[sample_idx*3+1]>llHapPairs[sample_idx*3+0]+.1 || llHapPairs[sample_idx*3+2]>llHapPairs[sample_idx*3+0]+.10) std::cerr << " *** ";
                std::cerr << std::endl;
           }
       
       
       }

       // estimate posterior
       const double lc = log(2.0);

       Q1[0] = llHapPairs[0*3+0];
       Q1[1] = lc+llHapPairs[0*3+1];
       Q1[2] = llHapPairs[0*3+2];
       
       for (int s=2;s<=numSamples;s++)
       {
       int sidx = s-1;

           Q2[0] = Q1[0]+llHapPairs[sidx*3+0];
           Q2[1] = addLogs(Q1[1]+llHapPairs[sidx*3+0],lc+Q1[0]+llHapPairs[sidx*3+1]);
           for (int k=2;k<=2*s-2;k++)
           {
                Q2[k] = addLogs(addLogs(Q1[k]+llHapPairs[sidx*3+0],lc+Q1[k-1]+llHapPairs[sidx*3+1]),Q1[k-2]+llHapPairs[sidx*3+2]);
           }
       Q2[2*s-1] = addLogs(lc+Q1[2*s-1-1]+llHapPairs[sidx*3+1],Q1[2*s-1-2]+llHapPairs[sidx*3+2]);
       Q2[2*s] = Q1[2*s-2]+llHapPairs[sidx*3+2];
           Q1.swap(Q2);

           if (DEBUG_CALLINDEL)
           {
               if (s==numSamples)
               {
                    std::cerr << "candidate haplotype " << *cand_it << " s: " << s;
                    double max_nonref = -HUGE_VAL;
                    int kidx=0;
                    for (int k=0;k<=2*s;k++)
                    {
                            std::cerr << " " << Q1[k];
                            if (k>0 && Q1[k]>max_nonref)
                            {
                                    kidx = k;
                                    max_nonref = Q1[k];
                            }
                    }
                    std::cerr << std::endl;
                    std::cerr << "max_k: " << kidx << " max_nonref: " << max_nonref << " diff: " << max_nonref - Q1[0] << " prior: " << prior[kidx] << std::endl;
               }
           }
       

       } 

       // Q1 holds the result

       double probFreq=-HUGE_VAL;
       double mapFreq = 0.0;
       double norm = -HUGE_VAL;
       for (int k=0;k<=2*numSamples;k++)
       {
           Q1[k] += prior[k];
           double post = Q1[k];
           norm = addLogs(post, norm);
           if (post>probFreq)
           {
               probFreq = post;
               mapFreq = double(k)/double(2*numSamples);
           }
       }
       for (int k=0;k<=2*numSamples;k++) 
       {
       Q1[k]-=norm;
           if (DEBUG_CALLINDEL) std::cerr << "post[" << k << "]: " << Q1[k] << std::endl;

       }


       double qual = -(Q1[0])*10.0/2.3056;
       qual_sample_lik *= (10.0/2.3056);
       if (qual >=10.0 || (qual<10.0 && qual>qual_sample_lik))
       {
           result[*cand_it] = CoalescentResult(qual, mapFreq, false, qual_sample_lik);
       } else 
       {
           // in this case it is based on the single sample caller.
           result[*cand_it] = CoalescentResult(qual_sample_lik, mapFreq, true, qual_sample_lik);
       }

       if (DEBUG_CALLINDEL) std::cerr << "Candidate haplotype " << *cand_it << " qual: " << qual << " freq: " << mapFreq << std::endl;
   }
   return result;
}



void DindelRealignWindow::printReadAlignments(int readIdx, std::ostream & out, int offset = 10, bool supportAlt = false)
{

    // print alignments against each candidate haplotype

   const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::printReadAlignments incomplete haplotype alignments.");
   DindelRead & read = (*m_pDindelReads)[readIdx];
   

   out << "ALIGNMENTS for read " << read.getID() << " sample: " << read.getSampleName() << " seq: " << read.getSequence() << " ref lik: " << hapReadAlignments[0][readIdx].logLik <<  std::endl;

   /*
   for (int h=0;h<int(haplotypes.size());h++)
   {

       const DindelHaplotype & haplotype = haplotypes[h];
       
       std::cerr << "\n HASH Haplotype " << h << std::endl;
       
       std::cerr << "\nreadkeys:" << std::endl;
       const std::vector<unsigned int> & readKeys = read.getHashKeys();
       for (int x=0;x<int(readKeys.size());x++) { 
        std::cerr << "\treadpos: " << x << " " << readKeys[x] << " hapHash lookup: " << haplotype.getHash().lookup(readKeys[x]) << std::endl;
       }
   }
   */

   int hstart = 0;
   if (supportAlt) hstart = 0;
   bool printRef = false;
   for (int h=int(haplotypes.size())-1;h>=hstart;h--)
   {
       if (supportAlt && h!=0 && (hapReadAlignments[h][readIdx].logLik<-10.0 || (h!=0 && hapReadAlignments[h][readIdx].logLik <= hapReadAlignments[0][readIdx].logLik))) continue;
       if (h!=0) printRef = true;
       if (h==0 && !printRef) break;
       const DindelHaplotype & haplotype = haplotypes[h];
       std::string pad(offset, ' ');


       // print variant annotations in the haplotype
       bool hasVars = false;
       out << pad;
       for (int b=0;b<haplotype.length();b++)
       {
           int rb = haplotype.getRefBase(b);

           if (rb<0) 
       {
           if (rb == SNP || rb == MULTINUCLEOTIDE_RUN) out << "S";
           else if (rb == INSERTION) out << "I";
           else out << -rb; 
           hasVars = true;
       }
       else 
           {
               if (b==0) out << " ";
               else
               {
                   // check for deletion
                   int prb = haplotype.getRefBase(b-1);
                   if (prb>=0 && rb>=0 && rb-prb>1) 
           {
               out << "D"; 
               hasVars = true;
           }
           else out << " ";
               }
           }
       }
       out << std::endl;

       // print haplotype
       std::ostringstream os;
       os << h;
       out << os.str() << std::string(offset-os.str().length(),' ') << haplotype.getSequence() << std::endl;


       // print read alignment
       const ReadHaplotypeAlignment & rha = hapReadAlignments[h][readIdx];

       int end = rha.hapPosLastReadBase;
       int start = end-read.length()+1;
       std::string readseq = read.getSequence();

       std::string rpad = "";
       if (start>0)
       {
          rpad = std::string(start, ' ');
       }

       out << pad << rpad;
       for (int b=start,r=0;b<=end && r<read.length();b++,r++)
       {
           if (b>=0)
           {
               if (b<haplotype.length())
               {
                   // on haplotype
                   if (haplotype.getSequence()[b]==readseq[r]) out << "."; else out << readseq[r];
               } 
           }
       }

       // output read haplotype alignment statistics. Still on same line

       out << " logLik: " << rha.logLik << " postProbLastReadBase: " << rha.postProbLastReadBase << " newpos: " << haplotype.getRefStart()+start << " start: " << start << " end: " << end;

       out << std::endl;
       out << std::endl;
   }
}



DindelRealignWindowResult DindelRealignWindow::estimateHaplotypeFrequencies(double minLogLikAlignToRef, double minLogLikAlignToAlt, bool capUsingMappingQuality, bool print = false)
{
   if (DINDEL_DEBUG) 
   {
       std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies START" << std::endl;
   }
    
   const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::estimateHaplotypeFrequencies incomplete haplotype alignments.");

   // init result
   DindelRealignWindowResult result(*this);


   // construct read X haplotype vector
   int numReads = (int) m_pDindelReads->size();
   const std::vector<DindelRead> & reads = (*m_pDindelReads);
   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik( numReads, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));
   for (int r=0;r<numReads;r++)
       for (int h=0;h<numHaps;h++)
       {
           if (haplotypes[h].isReference()) hrLik[r][h]=addLogs(hapReadAlignments[h][r].logLik, minLogLikAlignToRef);
            else hrLik[r][h]=addLogs(hapReadAlignments[h][r].logLik, minLogLikAlignToAlt);

           if (capUsingMappingQuality) hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));
        // if ((*m_pDindelReads)[r].getID()=="ind.0.1_473521_473690_d17") std::cerr << (*m_pDindelReads)[r].getID() << " r " << r << " h " << h << " lik: " << hrLik[r][h] << std::endl;
       }

   if (DEBUG_CALLINDEL)
   {
       for (int r=0;r<numReads;r++)
       {
            printReadAlignments(r, std::cerr,10,true);
       }
   }

   // now add haplotypes until we have obtained a maximal set of hapltoypes

   std::vector<int> added(haplotypes.size(),0);

   assert(haplotypes[0].isReference());
   added[0] = 1; // we always 'call' the reference haplotype


   size_t numAdded = 1;
   bool done = false;

   std::vector<double> addLL(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> coalQual(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> estFreqs(numHaps, 1e-10); // estimated haplotype frequencies
   std::vector< std::set<int> > addReads(numHaps*2); // first numHaps are forward, last are reverse

   estFreqs[0] = 1.0;

   typedef std::map<DindelVariant, HashMap<int, int> > VariantToHaplotype;
   VariantToHaplotype variantToHaplotype;

   for (int i=0;i<numHaps;i++) 
   {
       const std::vector<DindelVariant> & vars = haplotypes[i].getVariants();
       for (size_t x=0;x<vars.size();x++) variantToHaplotype[vars[x]][i]=1;
   }

   HapIdxToCoalescentResult hapCoalResults;
    
   int iteration = 0;
   while (!done)
   {
       if (numAdded == haplotypes.size()) break; // always keep at least one haplotype

       // compute for each haplotype the penalty for removing it.

       // only itialize kept haplotypes. 
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           addLL[h]=0.0;
           addReads[h].clear();
           addReads[h+numHaps].clear();
       }

       // setup haplotypes
       std::list<int> calledHaplotypes, candidateHaplotypes;
       std::vector<double> calledFreqs;

       for (int h=0;h<numHaps;h++)
       {
       // need the 
           if (added[h]) {
               calledHaplotypes.push_back(h);
               calledFreqs.push_back(estFreqs[h]);
           } else
           {

               // haplotype was not added. If in a previous iteration the coalQual was higher than 20, we will test it again to see if it explains more reads.
               if (iteration==0 || (iteration>=1 && coalQual[h]>=20.0)) candidateHaplotypes.push_back(h);
           }
       }
       // only test if there are still candidate haplotypes
       if (candidateHaplotypes.size()==0) break;
       
       hapCoalResults = computeHaplotypesCoalescent(calledHaplotypes, calledFreqs, candidateHaplotypes, minLogLikAlignToAlt);

    /*
    for (HapIdxToCoalescentResult::const_iterator it = hapCoalResults.begin();it!=hapCoalResults.end();it++) 
    {
        std::cerr << "iteration[" << iteration << "]: haplotype[" << it->first << "]: ";
        const std::vector<DindelVariant> & vars = haplotypes[it->first].getVariants();
        for (size_t v=0;v<vars.size();v++) std::cerr << " " << vars[v].getID();
        std::cerr << " qual: " << it->second.qual << " freq: " << it->second.mapFreq << std::endl;
    }
    */

       for (int r=0;r<numReads;r++)
       {
       //std::cerr << "READ " << r << std::endl;
           std::multiset<HL> hls;
           for (int h=0;h<numHaps;h++) if (added[h]) hls.insert(HL(hrLik[r][h],h));
           std::multiset<HL>::const_reverse_iterator best = hls.rbegin();

           double max_ll = -1000.0;
           if (best != hls.rend() && best->ll>max_ll) max_ll = best->ll;

       //bool print = false;
       for (int h=0;h<numHaps;h++) if (!added[h])
           {
               //double diff = hrLik[r][h]-max_ll;
           double min_ll_diff = HUGE_VAL;
           for (int h1=0;h1<numHaps;h1++) if (added[h1]) {
              double tdiff = hrLik[r][h]-hrLik[r][h1];
              if (tdiff<min_ll_diff) {
                  min_ll_diff = tdiff;
              }
           }
               //double diff = hrLik[r][h]-hrLik[r][0]; // only compare to reference haplotype?
           double diff = min_ll_diff;
               if (diff>0.0)
               {
                   addLL[h] += diff;
                   if (diff>2.0) 
           {
               if ((*m_pDindelReads)[r].isForward()) addReads[h].insert(r); else addReads[h+numHaps].insert(r);
           }
                   //if (false && realignParameters.showCallReads && h>0 && diff>2.0) print = true;
               }

           }
       
       }

       std::multiset<HL> addCandidates;
       std::map<int, std::multiset<HL> > orderedAddCandidates; // ordered by number of new variants added

        
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           // if we would remove haplotype h, how many variants would we loose?
       // determine current set of sequence variants


           // see which variants are unique to the haplotype
           
           int numUniqueVars = 0;
           double logPrior = 0.0;
       const std::vector<DindelVariant> & vars = haplotypes[h].getVariants();
           for (size_t x=0;x<vars.size();x++)
       {
          VariantToHaplotype::const_iterator iter = variantToHaplotype.find(vars[x]);
          //assert(iter != variantToHaplotype.end());
          if (! (iter != variantToHaplotype.end())) throw std::string("iter != variantToHaplotype.end()");
          bool varIsInOtherHaplotype=false;
          for (HashMap<int, int>::const_iterator ith = iter->second.begin();ith!=iter->second.end();ith++)
          {
              if (added[ith->first] && ith->first !=h)
              {
                  varIsInOtherHaplotype = true;
                  break;
              }
          }
          if (!varIsInOtherHaplotype) 
                  {
                      numUniqueVars++;
                      logPrior += log(vars[x].getPriorProb());
                  }
       }

           double diffLL = logPrior + addLL[h];
       if (DINDEL_DEBUG)
       {
        std::cerr << "estimate haplotype " << h << " numUniqueVars: " << numUniqueVars << " logPrior  " << logPrior << " addLL[h]: " << addLL[h] << std::endl;

       }
           
           addLL[h] = diffLL;
           //assert(hapCoalResults.find(h)!=hapCoalResults.end());
       double qual = coalQual[h];
       if (hapCoalResults.find(h)!=hapCoalResults.end())
       {
            qual = hapCoalResults[h].qual; // SET COALESCENT QUALITY
       }
           coalQual[h] = qual;
           orderedAddCandidates[numUniqueVars].insert(HL(qual, h));
       }
       
       // add the haplotype that adds the lowest number of variants with the highest LL change. 
       // this avoids the hitchhiking effect of one variant with only one supporting read when there is another haplotype that also improves the
       // score but increases the number of variants only by one.

       /*
       if (realignParameters.showCallReads)
       {
           for (size_t h=0;h<haplotypes.size();h++)
           {
               
               std::cerr << "Haplotype[" << h << "]: ";
               const std::vector<DindelVariant> & vars = haplotypes[h].getVariants();
               for (size_t x=0;x<vars.size();x++) std::cerr << " " << vars[x].getID();
               std::cerr << std::endl;
           }
       }
       */


       std::vector<DindelVariant> addedVariants;

       for (std::map<int, std::multiset<HL> >::const_iterator it1 = orderedAddCandidates.begin(); it1!=orderedAddCandidates.end();it1++)
       {
           std::multiset<HL>::const_reverse_iterator iter = it1->second.rbegin();
           if (iter != it1->second.rend())
           {
               done = true;
               if (iter->ll>=5.0) // NOTE HERE IS THE DECISION TO ADD A HAPLOTYPE BASED ON COALESCENT RESULTS
               {
                    int hapIdx = iter->idx;
                    
                    done = false;
                    added[hapIdx]=1;
                    numAdded++;
                estFreqs[hapIdx] = hapCoalResults[hapIdx].mapFreq;
            double fs = 0.0;
            for (int y=0;y<numHaps;y++) fs += estFreqs[y];
            for (int y=0;y<numHaps;y++) estFreqs[y]/=fs;

                    bool indelAdded = false;

                    const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants();
                    for (size_t x=0;x<vars.size();x++) if (result.variantInference.find(vars[x])==result.variantInference.end())
                    {
            
                        if (vars[x].getType()=="INDEL") indelAdded = true;
                        // this is a variant that was not called yet from other haplotypes
                        std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(vars[x],DindelRealignWindowResult::Inference()));
                        DindelRealignWindowResult::Inference & varInf = it_pair.first->second;

                        varInf.numReadsForward=int(addReads[hapIdx].size());
                        varInf.numReadsReverse=int(addReads[hapIdx+numHaps].size());
                        varInf.numReadsForwardZeroMismatch=0;
                        varInf.numReadsReverseZeroMismatch=0;
            varInf.numUnmapped = 0;
                        varInf.numLibraries = 0;
                        varInf.numReadNames = 0;



                        varInf.qual=iter->ll; //addLL[hapIdx]/.23026;
                        varInf.freq=hapCoalResults[hapIdx].mapFreq;
                        varInf.isSingleSampleCall = hapCoalResults[hapIdx].singleSampleCall;
            varInf.singleSampleQual = hapCoalResults[hapIdx].singleSampleQual;
                                    // find out which reads mapped without mismatches..
                        HashMap<std::string, int> libraries, readnames;
                        for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r!= addReads[hapIdx].end(); r++)
            {
                if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsForwardZeroMismatch+=1;
                if ((*m_pDindelReads)[*r].isUnmapped()) varInf.numUnmapped++;
                                if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped)  varInf.addDistanceToHistogram(haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()));
 
                varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][*r].logLik/double(reads[*r].length()));
                //if (indelAdded) std::cerr << "indel added dist: " << haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()) << std::endl;
                                varInf.addMapQToHistogram(reads[*r].getMappingQual());
                                libraries[reads[*r].getLibraryName()]=1;
                                readnames[reads[*r].getID()]=1;
                if (DINDEL_DEBUG) {
                    std::cerr << "==================== CHECK DIST HISTO HAPLOTYPE " << hapIdx << " with INDEL " << std::endl;
                    std::cerr << "  variants: ";
                    std::cerr << "NMM:  " << hapReadAlignments[hapIdx][*r].nmm << "isUngapped: " << hapReadAlignments[hapIdx][*r].isUngapped << " >minLogLikAligntoAlt: " << int(hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt) << "ll: " << hapReadAlignments[hapIdx][*r].logLik << " minLogLikAlignToAlt: " << minLogLikAlignToAlt << std::endl;

                    for (size_t xx=0;xx<vars.size();xx++) std::cerr << " " << vars[xx].getID(); std::cerr << std::endl;

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);
                    if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped) std::cerr << "HISTDIST ADD FOR READ "<< *r << " == " << (haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length())) << std::endl;
     
                    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

                }
                }
                        for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r!= addReads[hapIdx+numHaps].end(); r++) 
            {
                if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsReverseZeroMismatch+=1;
                if ((*m_pDindelReads)[*r].isUnmapped()) varInf.numUnmapped++;
                                if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped)  varInf.addDistanceToHistogram(haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()));
                                varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][*r].logLik/double(reads[*r].length()));
                                varInf.addMapQToHistogram(reads[*r].getMappingQual());
                // if (indelAdded) std::cerr << "indel added dist: " << haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()) << std::endl;
                                libraries[reads[*r].getLibraryName()]=1;
                                readnames[reads[*r].getID()]=1;
                if (DINDEL_DEBUG) {
                    std::cerr << "==================== CHECK DIST HISTO HAPLOTYPE " << hapIdx << " with INDEL " << std::endl;
                    std::cerr << "  variants: ";
                    std::cerr << "NMM:  " << hapReadAlignments[hapIdx][*r].nmm << "isUngapped: " << hapReadAlignments[hapIdx][*r].isUngapped << " >minLogLikAligntoAlt: " << int(hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt) << "ll: " << hapReadAlignments[hapIdx][*r].logLik << " minLogLikAlignToAlt: " << minLogLikAlignToAlt << std::endl;


                    for (size_t xx=0;xx<vars.size();xx++) std::cerr << " " << vars[xx].getID(); std::cerr << std::endl;

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);
                        if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped) std::cerr << "HISTDIST ADD FOR READ "<< *r << " == " << (haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length())) << std::endl; 
                    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

                }

                        }
                        varInf.strandBias = DindelRealignWindowResult::Inference::computeStrandBias(varInf.numReadsForward, varInf.numReadsReverse);
                        varInf.numLibraries = int(libraries.size());
                        varInf.numReadNames = int(readnames.size());
            }

                    if (indelAdded && realignParameters.showCallReads!=0)
                    {
                        std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies Adding haplotype " << hapIdx << " with INDEL " << std::endl;
                        std::cerr << "  variants: ";
                        for (size_t x=0;x<vars.size();x++) std::cerr << " " << vars[x].getID(); std::cerr << std::endl;

                        for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r != addReads[hapIdx].end();r++)
                            printReadAlignments(*r, std::cerr);

                        for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r != addReads[hapIdx+numHaps].end();r++)
                            printReadAlignments(*r, std::cerr);



                    }

                    // only add one haplotype at a time
                    break;
               }
           }
    }
       iteration++;
   }




   if (!QUIET) std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies Added " << numAdded << " out of " << haplotypes.size() << " haplotypes." << std::endl;
   if (realignParameters.showCallReads && print)
   {
       for (int h=0;h<numHaps;h++)
       {
               std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies final penalty haplotype[" <<h << "]: " << addLL[h] << std::endl;
       }

   }


   // see which variants were not called
   std::vector<DindelVariant> notCalled;
   
   for (VariantToHaplotype::const_iterator iter=variantToHaplotype.begin();iter!=variantToHaplotype.end();iter++)
   {
       if (result.variantInference.find(iter->first)==result.variantInference.end())
       {
            // not called
            double maxLL = -HUGE_VAL;
            int maxIdx=-1;
            int hapIdx;
            for (HashMap<int, int>::const_iterator ith = iter->second.begin();ith!=iter->second.end();ith++)
            {
               hapIdx = ith->first;
               if (addLL[hapIdx]>maxLL)
               {
                   maxIdx = hapIdx;
                   maxLL = addLL[hapIdx];
               }
            }
            assert(maxIdx!=-1);
            hapIdx = maxIdx;

            std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(iter->first,DindelRealignWindowResult::Inference()));
            DindelRealignWindowResult::Inference & varInf = it_pair.first->second;

            varInf.numReadsForward=int(addReads[hapIdx].size());
            varInf.numReadsReverse=int(addReads[hapIdx+numHaps].size());
            varInf.numReadsForwardZeroMismatch=0;
            varInf.numReadsReverseZeroMismatch=0;

            varInf.qual=hapCoalResults[hapIdx].qual; //addLL[hapIdx]/.23026;
            varInf.freq=-1.0;
            // find out which reads mapped without mismatches..
            for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r!= addReads[hapIdx].end(); r++) if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsForwardZeroMismatch+=1;
            for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r!= addReads[hapIdx+numHaps].end(); r++) if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsReverseZeroMismatch+=1;


       }

   }
   
  // assign result


   result.haplotypes = haplotypes;
   result.haplotypeFrequencies = std::vector<double>(numHaps, 0.0);
   size_t totReads = 0; for (size_t x=0;x<addReads.size();x++) totReads += addReads[x].size();
   
   for (int h=0;h<numHaps;h++)
   {
       result.haplotypeFrequencies[h] = estFreqs[h];
   }

   result.numHapSpecificReads = (int) totReads;
   result.numHaplotypes = (int) numAdded-1; // -1 makes sure we don't count the reference haplotype. CHANGE when we are doing proper genotyping, because then we align each read to each haplotype
   return result;
   
}


void DindelRealignWindow::doEM(const std::vector< std::vector<double> > & hrLik, const std::vector<int> & calledHaplotypes, std::vector<double> & haplotypeFrequencies)
{




    // first remove uncalled haplotypes
    size_t nh=0;
    for(size_t h=0;h < calledHaplotypes.size(); ++h) if (calledHaplotypes[h]) nh++;
    size_t nr = hrLik.size();


    std::vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods1
    
    size_t hidx = 0;
    for (size_t h=0;h<calledHaplotypes.size();h++) if (calledHaplotypes[h]) {
        
        for (size_t r=0;r<nr;r++) {
            rl[hidx*nr+r]=hrLik[r][h];
        }
        ++hidx;
    }

    if (DINDEL_DEBUG) std::cerr << "doEM: nh: " << nh << " nr: " << nr << std::endl;
    

    std::vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
    std::vector<double> pi(nh); // log of haplotype frequencies
    std::vector<double> nk(nh,0.0); // counts for each haplotype

    std::vector<double> hapFreqs=nk;

    // initialize frequencies
    for (size_t h=0;h<nh;h++) pi[h]=1.0/double(nh);

    // initialize expectations of indicator variables
    for (size_t h=0;h<nh;h++) for (size_t r=0;r<nr;r++) {
            z[h*nr+r]=0.5;
    }


    bool converged=false;
    double tol=this->realignParameters.EMtol;

    double eOld=-HUGE_VAL, eNew;

    int iter=0;
    while (!converged) {

            // compute expectation of indicator variables
            for (size_t h=0;h<nh;h++) nk[h]=0.0;

            int idx=0;
            for (size_t r=0;r<nr;r++) {
                    double lognorm=-HUGE_VAL;
                    // compute responsibilities
                    for (size_t h=0;h<nh;h++) {
                            z[h*nr+r]=pi[h]+(rl[h*nr+r]);
                            lognorm=addLogs(lognorm, z[h*nr+r]);
                    }
                    // normalize and compute counts
                    for (size_t h=0;h<nh;h++) {
                            z[nr*h+r]-=lognorm;
                            z[nr*h+r]=exp(z[nr*h+r]);

                            nk[h]+=z[nr*h+r];
                    }
            }

            // compute frequencies

            for (size_t h=0;h<nh;h++) {
                    double nph=nk[h]/nr;
                    pi[h]=log(nph);
            }


            idx=0;
            eNew=0.0;
            for (size_t h=0;h<nh;h++)
            {
                for (size_t r=0;r<nr;r++)
                {
                    // compute responsibilities
                    eNew+=z[idx]*( pi[h]+rl[idx]);
                    idx++;
                }
            }
            if (DINDEL_DEBUG) std::cerr << " EM iter: " << iter << " " << eNew << " eOld-eNew: " <<  eOld-eNew << std::endl;
            //
            //assert (eOld-eNew<1e-10);
            converged=(fabs(eOld-eNew))<tol || iter>realignParameters.EMmaxiter;

            eOld=eNew;


            iter++;
    }

    haplotypeFrequencies = std::vector<double>(calledHaplotypes.size(),0.0);

    hidx = 0;
    for(size_t h=0;h < calledHaplotypes.size();h++)
    {
        if (calledHaplotypes[h]) haplotypeFrequencies[h] = exp(pi[hidx++]); else haplotypeFrequencies[h] = 0.0;
    }

    if(DINDEL_DEBUG)
    {
        std::cerr  << "DindelRealignWindow::doEM haplotype Frequencies: ";
        hidx = 0;
        for (size_t h=0;h<calledHaplotypes.size();h++) std::cerr << " [ " << h << " " << haplotypeFrequencies[h] << "];";
        std::cerr << std::endl;
    }
    
    
}

DindelRealignWindowResult DindelRealignWindow::estimateHaplotypeFrequenciesModelSelection(double minLogLikAlignToRef, double minLogLikAlignToAlt, bool capUsingMappingQuality, bool print = false)
{
   if (DINDEL_DEBUG)
   {
       std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies START" << std::endl;
   }

   const std::vector<DindelHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::estimateHaplotypeFrequencies incomplete haplotype alignments.");

   // init result
   DindelRealignWindowResult result(*this);


   // construct read X haplotype vector
   int numReads = (int) m_pDindelReads->size();
   const std::vector<DindelRead> & reads = (*m_pDindelReads);
   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik( numReads, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));
   for (int r=0;r<numReads;r++)
       for (int h=0;h<numHaps;h++)
       {
           if (haplotypes[h].isReference()) hrLik[r][h]=addLogs(hapReadAlignments[h][r].logLik, minLogLikAlignToRef);
            else hrLik[r][h]=addLogs(hapReadAlignments[h][r].logLik, minLogLikAlignToAlt);

           if (capUsingMappingQuality) hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));
        // if ((*m_pDindelReads)[r].getID()=="ind.0.1_473521_473690_d17") std::cerr << (*m_pDindelReads)[r].getID() << " r " << r << " h " << h << " lik: " << hrLik[r][h] << std::endl;
       }

   if (DEBUG_CALLINDEL)
   {
       for (int r=0;r<numReads;r++)
       {
            printReadAlignments(r, std::cerr,10,true);
       }
   }

   // now add haplotypes until we have obtained a maximal set of hapltoypes

   std::vector<int> added(haplotypes.size(),0);

   assert(haplotypes[0].isReference());
   added[0] = 1; // we always 'call' the reference haplotype


   size_t numAdded = 1;
   bool done = false;

   std::vector<double> addLL(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> posteriorAddLL(numHaps, 0.0); // loglikelihood gained by adding haplotypes but penalized for the number of variants that will be added.
   std::vector<double> haplotypeQual(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> haplotypeFrequencies(numHaps, 0.0);
   
   std::vector< std::set<int> > addReads(numHaps*2); // first numHaps are forward, last are reverse

   // map variants to haplotypes
   typedef std::map<DindelVariant, HashMap<int, int> > VariantToHaplotype;
   VariantToHaplotype variantToHaplotype;

   for (int i=0;i<numHaps;i++)
   {
       const std::vector<DindelVariant> & vars = haplotypes[i].getVariants();
       for (size_t x=0;x<vars.size();x++) variantToHaplotype[vars[x]][i]=1;
   }

   
   int iteration = 0;
   while (!done)
   {
       if (numAdded == haplotypes.size()) break; // always keep at least one haplotype

       // compute for each haplotype the penalty for removing it.

       // only itialize kept haplotypes.
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           addLL[h]=0.0;
               posteriorAddLL[h]=-HUGE_VAL;
           addReads[h].clear();
           addReads[h+numHaps].clear();
       }

       // setup haplotypes
       std::list<int> calledHaplotypes, candidateHaplotypes;
      

       for (int h=0;h<numHaps;h++)
       {
       // need the
           if (added[h]) {
               calledHaplotypes.push_back(h);
           } else
           {

               // haplotype was not added. If in a previous iteration the coalQual was higher than 20, we will test it again to see if it explains more reads.
               // CHECK IF THIS STILL MAKES SENSE FIXME
               if (iteration==0 || (iteration>=1 && haplotypeQual[h]>=20.0)) candidateHaplotypes.push_back(h);
           }
       }
       // only test if there are still candidate haplotypes
       if (candidateHaplotypes.size()==0) break;


     
       for (int r=0;r<numReads;r++)
       {
       //std::cerr << "READ " << r << std::endl;
           std::multiset<HL> hls;
           for (int h=0;h<numHaps;h++) if (added[h]) hls.insert(HL(hrLik[r][h],h));
           std::multiset<HL>::const_reverse_iterator best = hls.rbegin();

           double max_ll = -1000.0;
           if (best != hls.rend() && best->ll>max_ll) max_ll = best->ll;

       //bool print = false;
       for (int h=0;h<numHaps;h++) if (!added[h])
           {
               //double diff = hrLik[r][h]-max_ll;

               // find the called haplotype that gives the best likelihood.
           double min_ll_diff = HUGE_VAL;
           for (int h1=0;h1<numHaps;h1++) if (added[h1]) {
              double tdiff = hrLik[r][h]-hrLik[r][h1];
              if (tdiff<min_ll_diff) {
                  min_ll_diff = tdiff;
              }
           }
               //double diff = hrLik[r][h]-hrLik[r][0]; // only compare to reference haplotype?
           double diff = min_ll_diff;
               if (diff>0.0)
               {
                   addLL[h] += diff;
                   if (diff>2.0)
           {
               if ((*m_pDindelReads)[r].isForward()) addReads[h].insert(r); else addReads[h+numHaps].insert(r);
           }
                   //if (false && realignParameters.showCallReads && h>0 && diff>2.0) print = true;
               }

           }

       }
       std::map<int, std::multiset<HL> > orderedAddCandidates; // ordered by number of new variants added


       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           // if we would add haplotype h, what would be the increase in the likelihood?
       // determine current set of sequence variants


           // see which variants are unique to the haplotype

           int numUniqueVars = 0;
           double logPrior = 0.0;
       const std::vector<DindelVariant> & vars = haplotypes[h].getVariants();
           for (size_t x=0;x<vars.size();x++)
       {
          VariantToHaplotype::const_iterator iter = variantToHaplotype.find(vars[x]);
          //assert(iter != variantToHaplotype.end());
          if (! (iter != variantToHaplotype.end())) throw std::string("iter != variantToHaplotype.end()");
          bool varIsInOtherHaplotype=false;
          for (HashMap<int, int>::const_iterator ith = iter->second.begin();ith!=iter->second.end();ith++)
          {
              if (added[ith->first] && ith->first !=h)
              {
                  varIsInOtherHaplotype = true;
                  break;
              }
          }
          if (!varIsInOtherHaplotype)
                  {
                      numUniqueVars++;
                      logPrior += log(vars[x].getPriorProb());
                  }
       }

           double qual = logPrior + addLL[h];
       if (DINDEL_DEBUG)
       {
        std::cerr << "estimate haplotype " << h << " numUniqueVars: " << numUniqueVars << " logPrior  " << logPrior << " haplotype qual: " << qual << std::endl;

       }

           posteriorAddLL[h] = qual;
           haplotypeQual[h] = qual;

           // NOTE/FIXME, prior should take care of hitchhiking effects. CHECK IF THIS IS TRUE....
           numUniqueVars = 1;
           orderedAddCandidates[numUniqueVars].insert(HL(qual, h));
       }

       // add the haplotype that adds the lowest number of variants with the highest LL change.
       // this avoids the hitchhiking effect of one variant with only one supporting read when there is another haplotype that also improves the
       // score but increases the number of variants only by one.
       // SHOULD BE COMPENSATED FOR BY logPRIOR? WHAT IF VARIANTS ARE LINKED?

       /*
       if (realignParameters.showCallReads)
       {
           for (size_t h=0;h<haplotypes.size();h++)
           {

               std::cerr << "Haplotype[" << h << "]: ";
               const std::vector<DindelVariant> & vars = haplotypes[h].getVariants();
               for (size_t x=0;x<vars.size();x++) std::cerr << " " << vars[x].getID();
               std::cerr << std::endl;
           }
       }
       */


  
       for (std::map<int, std::multiset<HL> >::const_iterator it1 = orderedAddCandidates.begin(); it1!=orderedAddCandidates.end();it1++)
       {
           std::multiset<HL>::const_reverse_iterator iter = it1->second.rbegin();
           if (iter != it1->second.rend())
           {
               done = true;
               if (iter->ll>=5.0) // NOTE HERE IS THE DECISION TO ADD A HAPLOTYPE BASED ON COALESCENT RESULTS
               {
                    int hapIdx = iter->idx;

                    done = false;
                    added[hapIdx]=1;
                    numAdded++;
             
                    bool indelAdded = false;

                    const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants();
                    for (size_t x=0;x<vars.size();x++) if (result.variantInference.find(vars[x])==result.variantInference.end())
                    {

                        if (vars[x].getType()=="INDEL") indelAdded = true;
                        // this is a variant that was not called yet from other haplotypes
                        std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(vars[x],DindelRealignWindowResult::Inference()));
                        DindelRealignWindowResult::Inference & varInf = it_pair.first->second;

                        varInf.numReadsForward=int(addReads[hapIdx].size());
                        varInf.numReadsReverse=int(addReads[hapIdx+numHaps].size());
                        varInf.numReadsForwardZeroMismatch=0;
                        varInf.numReadsReverseZeroMismatch=0;
            varInf.numUnmapped = 0;
                        varInf.numLibraries = 0;
                        varInf.numReadNames = 0;



                        varInf.qual=iter->ll; //addLL[hapIdx]/.23026;
                        varInf.freq=0.0;
                        varInf.isSingleSampleCall = 0;
            varInf.singleSampleQual = -1;
                                    // find out which reads mapped without mismatches..
                        HashMap<std::string, int> libraries, readnames;
                        for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r!= addReads[hapIdx].end(); r++)
            {
                if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsForwardZeroMismatch+=1;
                if ((*m_pDindelReads)[*r].isUnmapped()) varInf.numUnmapped++;
                                if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped)  varInf.addDistanceToHistogram(haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()));

                varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][*r].logLik/double(reads[*r].length()));
                //if (indelAdded) std::cerr << "indel added dist: " << haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()) << std::endl;
                                varInf.addMapQToHistogram(reads[*r].getMappingQual());
                                libraries[reads[*r].getLibraryName()]=1;
                                readnames[reads[*r].getID()]=1;
                if (DINDEL_DEBUG) {
                    std::cerr << "==================== CHECK DIST HISTO HAPLOTYPE " << hapIdx << " with INDEL " << std::endl;
                    std::cerr << "  variants: ";
                    std::cerr << "NMM:  " << hapReadAlignments[hapIdx][*r].nmm << "isUngapped: " << hapReadAlignments[hapIdx][*r].isUngapped << " >minLogLikAligntoAlt: " << int(hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt) << "ll: " << hapReadAlignments[hapIdx][*r].logLik << " minLogLikAlignToAlt: " << minLogLikAlignToAlt << std::endl;

                    for (size_t xx=0;xx<vars.size();xx++) std::cerr << " " << vars[xx].getID(); std::cerr << std::endl;

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);
                    if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped) std::cerr << "HISTDIST ADD FOR READ "<< *r << " == " << (haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length())) << std::endl;

                    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

                }
                }
                        for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r!= addReads[hapIdx+numHaps].end(); r++)
            {
                if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsReverseZeroMismatch+=1;
                if ((*m_pDindelReads)[*r].isUnmapped()) varInf.numUnmapped++;
                                if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped)  varInf.addDistanceToHistogram(haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()));
                                varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][*r].logLik/double(reads[*r].length()));
                                varInf.addMapQToHistogram(reads[*r].getMappingQual());
                // if (indelAdded) std::cerr << "indel added dist: " << haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length()) << std::endl;
                                libraries[reads[*r].getLibraryName()]=1;
                                readnames[reads[*r].getID()]=1;
                if (DINDEL_DEBUG) {
                    std::cerr << "==================== CHECK DIST HISTO HAPLOTYPE " << hapIdx << " with INDEL " << std::endl;
                    std::cerr << "  variants: ";
                    std::cerr << "NMM:  " << hapReadAlignments[hapIdx][*r].nmm << "isUngapped: " << hapReadAlignments[hapIdx][*r].isUngapped << " >minLogLikAligntoAlt: " << int(hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt) << "ll: " << hapReadAlignments[hapIdx][*r].logLik << " minLogLikAlignToAlt: " << minLogLikAlignToAlt << std::endl;


                    for (size_t xx=0;xx<vars.size();xx++) std::cerr << " " << vars[xx].getID(); std::cerr << std::endl;

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);

                        std::cerr << "FOR read " << *r << std::endl;
                        printReadAlignments(*r, std::cerr);
                        if (hapReadAlignments[hapIdx][*r].logLik>minLogLikAlignToAlt && hapReadAlignments[hapIdx][*r].nmm<=2 && hapReadAlignments[hapIdx][*r].nmm>=0 && hapReadAlignments[hapIdx][*r].isUngapped) std::cerr << "HISTDIST ADD FOR READ "<< *r << " == " << (haplotypes[hapIdx].getClosestDistance(vars[x],hapReadAlignments[hapIdx][*r].hapPosLastReadBase, hapReadAlignments[hapIdx][*r].hapPosLastReadBase-reads[*r].length())) << std::endl;
                    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

                }

                        }
                        varInf.strandBias = DindelRealignWindowResult::Inference::computeStrandBias(varInf.numReadsForward, varInf.numReadsReverse);
                        varInf.numLibraries = int(libraries.size());
                        varInf.numReadNames = int(readnames.size());
            }

                    if (indelAdded && realignParameters.showCallReads!=0)
                    {
                        std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies Adding haplotype " << hapIdx << " with INDEL " << std::endl;
                        std::cerr << "  variants: ";
                        for (size_t x=0;x<vars.size();x++) std::cerr << " " << vars[x].getID(); std::cerr << std::endl;

                        for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r != addReads[hapIdx].end();r++)
                            printReadAlignments(*r, std::cerr);

                        for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r != addReads[hapIdx+numHaps].end();r++)
                            printReadAlignments(*r, std::cerr);



                    }

                    // only add one haplotype at a time
                    break;
               }
           }
    }
       iteration++;
   }

   // estimate haplotype frequencies for called haplotypes using EM.

   doEM(hrLik, added, haplotypeFrequencies);
   

   // set variant frequencies
   for (VariantToHaplotype::const_iterator iter=variantToHaplotype.begin();iter!=variantToHaplotype.end();iter++)
   {
       
       DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(iter->first);
       if (vit!=result.variantInference.end())
       {
            for (HashMap<int, int>::const_iterator hit = iter->second.begin(); hit != iter->second.end(); hit++)
            {
                if (hit->second==1) vit->second.freq += haplotypeFrequencies[hit->first];
            }
       }
   }

   if (!QUIET) std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies Added " << numAdded << " out of " << haplotypes.size() << " haplotypes." << std::endl;

   if (realignParameters.showCallReads && print)
   {
       for (int h=0;h<numHaps;h++)
       {
               std::cerr << "DindelRealignWindow::estimateHaplotypeFrequencies final penalty haplotype[" <<h << "]: " << posteriorAddLL[h] << std::endl;
       }

   }


   // Uncalled variants
   for (VariantToHaplotype::const_iterator iter=variantToHaplotype.begin();iter!=variantToHaplotype.end();iter++)
   {
       if (result.variantInference.find(iter->first)==result.variantInference.end())
       {
            // not called
            double maxLL = -HUGE_VAL;
            int maxIdx=-1;
            int hapIdx;
            for (HashMap<int, int>::const_iterator ith = iter->second.begin();ith!=iter->second.end();ith++)
            {
               hapIdx = ith->first;
               if (addLL[hapIdx]>maxLL)
               {
                   maxIdx = hapIdx;
                   maxLL = addLL[hapIdx];
               }
            }
            assert(maxIdx!=-1);
            hapIdx = maxIdx;

            std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(iter->first,DindelRealignWindowResult::Inference()));
            DindelRealignWindowResult::Inference & varInf = it_pair.first->second;

            varInf.numReadsForward=int(addReads[hapIdx].size());
            varInf.numReadsReverse=int(addReads[hapIdx+numHaps].size());
            varInf.numReadsForwardZeroMismatch=0;
            varInf.numReadsReverseZeroMismatch=0;

            varInf.qual=haplotypeQual[hapIdx]; //addLL[hapIdx]/.23026;
            varInf.freq=0.0;
            // find out which reads mapped without mismatches..
            for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r!= addReads[hapIdx].end(); r++) if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsForwardZeroMismatch+=1;
            for (std::set<int>::const_iterator r = addReads[hapIdx+numHaps].begin(); r!= addReads[hapIdx+numHaps].end(); r++) if (hapReadAlignments[hapIdx][*r].nmm==0) varInf.numReadsReverseZeroMismatch+=1;


       }

   }

  // assign result


   result.haplotypes = haplotypes;
   result.haplotypeFrequencies = haplotypeFrequencies;
   size_t totReads = 0; for (size_t x=0;x<addReads.size();x++) totReads += addReads[x].size();

   

   result.numHapSpecificReads = (int) totReads;
   result.numHaplotypes = (int) numAdded-1; // -1 makes sure we don't count the reference haplotype. CHANGE when we are doing proper genotyping, because then we align each read to each haplotype
   return result;

}

WindowFile::WindowFile(const std::string& fileName)
{
    m_isOpen = false;
    m_fileName = fileName;
    m_fileHandle.open(fileName.c_str());
    if (!m_fileHandle.is_open())
    {
        std::string msg = "WindowFile::WindowFile::cannot_open_file_";
        msg.append(fileName);
        throw msg;
    }
    m_isOpen = true;
    m_index = 0;

    m_empty = std::vector<DindelVariant>();
}

std::vector<DindelVariant> WindowFile::getNextWindow()
{
    std::string line;
    std::getline(m_fileHandle, line);
    if (line.empty()) return m_empty;
    
    // update line counter
    m_index++;

    std::istringstream is(line);
  
    std::vector<DindelVariant> variants;

    while (!is.eof()) {
        std::string varstring;
        is >> varstring;
    
        DindelVariant var;
        bool success=DindelVariant::variantFromWindowString(varstring, var);
        if (!success) 
        {
            std::cerr << "WindowFile::getNextWindow Could not parse variant " + varstring + " in line " << m_index+1 << " of window file " << m_fileName << std::endl;
        }
        else
        {
            variants.push_back(var);
        }
    }

    return variants;
}

WindowFile::~WindowFile()
{
    if (m_isOpen)
    {
        m_fileHandle.close();
    }

}
/*
 *
 * REALIGNALGORITHM
 *
 */

DindelRealignParameters::DindelRealignParameters()
{
     initFromString(getDefaultParameters());
}

DindelRealignParameters::DindelRealignParameters(const std::string & paramString)
{
    initFromString(getDefaultParameters());
    // override with user parameters
    initFromString(paramString);
}

std::string DindelRealignParameters::getDefaultParameters() const
{
     std::string paramString = "genotyping:0,maxNumReads:100000,maxNumReadsWindow:100000,showCallReads:0,minNumHaplotypeOverlaps:0,maxNumCandidatesPerWindow:32,windowReadBuffer:500,minVariantSep:10,haplotypeWidth:60,minCandidateAlleleCount:0,probSNP:0.001,probINDEL:0.0001,probMNP:0.00001";
     paramString += ",maxMappingQuality:80,addSNPMaxSNPs:3,addSNPMaxMismatches:3,addSNPMinMappingQual:30,addSNPMinBaseQual:20,addSNPMinCount:2,minPostProbLastReadBaseForUngapped:0.95";
     paramString += ",singleSampleHetThreshold:20,singleSampleHomThreshold:20,EMtol:0.0001,EMmaxiter:200,doEM:1";
     return paramString; 
}      


void DindelRealignParameters::checkAndInit()
{
    bool print = false; 
    this->minLogLikAlignToAlt=-double(this->maxMappingQuality)*.2302585-log(double(DINDEL_HMM_BANDWIDTH)); // note last term accounts for base prior in HMM.
    this->minLogLikAlignToRef=this->minLogLikAlignToAlt;
    this->addSNPMinLogBaseQual=log(pow(10.0, -this->addSNPMinBaseQual/10.0));

    std::ostringstream os;

    os <<  "probSNP:" << variantPriors.m_probSNP;
    os << ",probINDEL:" << variantPriors.m_probINDEL;
    os << ",probMNP:" << variantPriors.m_probMNP;
    os << ",minVariantSep:" <<  minVariantSep;
    os << ",windowReadBuffer:" <<  windowReadBuffer;
    os << ",haplotypeWidth:" <<haplotypeWidth;
    os << ",minCandidateAlleleCount:" <<minCandidateAlleleCount;
    os << ",genotyping:" << genotyping;
    os << ",showCallReads:" << showCallReads;
    os << ",minNumHaplotypeOverlaps:" << minNumHaplotypeOverlaps;
    os << ",maxMappingQuality:" << maxMappingQuality;
    os << ",maxMappingQuality:" << maxMappingQuality;
    os << ",maxNumReads:" << maxNumReads;
    os << ",addSNPMinMappingQual:" << addSNPMinMappingQual;
    os << ",addSNPMinBaseQual:" << addSNPMinBaseQual;
    os << ",addSNPMinCount:" << addSNPMinCount;
    os << ",addSNPMaxSNPs:" << addSNPMaxSNPs;
    os << ",minPostProbLastReadBaseForUngapped:" << minPostProbLastReadBaseForUngapped;
    os << ",maxNumCandidatesPerWindow:" << maxNumCandidatesPerWindow;
    os << ",maxNumReadsWindow:" << maxNumReadsWindow;
    os << ",singleSampleHetThreshold:" << int(singleSampleHetThreshold);
    os << ",singleSampleHomThreshold:" << int(singleSampleHomThreshold);
    os << ",EMtol:" << EMtol;
    os << ",EMmaxiter:" << int(EMmaxiter);
    os << ",doEM:" << int(doEM);
    if (print) 
    {
        std::cout << "Realignment parameters:\n";
        this->m_paramString = os.str();
        for (size_t x=0;x<m_paramString.size();x++)
        {
        char c = m_paramString[x];
        if (c==':') std::cout << "\t";
        else if (c== ',') std::cout << "\n";
        else std::cout << c;
        }
        std::cout << std::endl;
    }

}

void DindelRealignParameters::initFromString(const std::string& paramString)
{
    assert(!paramString.empty());

    // check that it contains no spaces
    for (size_t x=0;x<paramString.size();x++) if (paramString[x]==' ' || paramString[x]=='\t')
    {
        std::cerr << "Parameter string cannot contain spaces or tabs." << std::endl;
        exit(1);
    }

    std::vector<std::string> pairs = SplitString(paramString,',');

    for (size_t p=0;p<pairs.size();p++)
    {
        std::vector<std::string> pair = SplitString(pairs[p],':');

        if (pair.size()==1)
        {
            // don't have any yet
        std::cerr << "Unrecognized flag" << std::endl;
        exit(1);
            
    }
        else if (pair.size()==2)
        {
            const std::string & k = pair[0];
            const std::string & v = pair[1];
            bool fail=false;
            double minProbVariant = 1e-10;
            double maxProbVariant = 1.0;

            if (k == "probSNP") { if (!from_string<double>(variantPriors.m_probSNP,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "probINDEL") { if (!from_string<double>(variantPriors.m_probINDEL,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "probMNP") { if (!from_string<double>(variantPriors.m_probMNP,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "minVariantSep") { if (!from_string<int>(minVariantSep,1, 100, v, std::dec)) fail = true; }
            else if (k == "windowReadBuffer") { if (!from_string<int>(windowReadBuffer, 10, 10000, v, std::dec)) fail = true; }
            else if (k == "haplotypeWidth") { if (!from_string<int>(haplotypeWidth, 20, 100, v, std::dec)) fail = true; }
        else if (k == "minCandidateAlleleCount") { if (!from_string<int>(minCandidateAlleleCount,0, 1000, v, std::dec)) fail = true; }
            else if (k == "maxMappingQuality") { if (!from_string<int>(maxMappingQuality,40,100,v, std::dec)) fail = true; }
            else if (k == "addSNPMaxMismatches") { if (!from_string<int>(addSNPMaxMismatches,1,10, v, std::dec)) fail = true; }
            else if (k == "addSNPMinMappingQual") { if (!from_string<double>(addSNPMinMappingQual,5.,50., v, std::dec)) fail = true; }
            else if (k == "addSNPMinBaseQual") { if (!from_string<int>(addSNPMinBaseQual,0, 20, v, std::dec)) fail = true; }
            else if (k == "addSNPMinCount") { if (!from_string<size_t>(addSNPMinCount,1, 20, v, std::dec)) fail = true; }
            else if (k == "addSNPMaxSNPs") { if (!from_string<int>(addSNPMaxSNPs,1, 6, v, std::dec)) fail = true; }
        else if (k == "minPostProbLastReadBaseForUngapped") { if (!from_string<double>(minPostProbLastReadBaseForUngapped,0.0, 1.0, v, std::dec)) fail = true; }
            else if (k == "genotyping") { if (!from_string<int>(genotyping,0, 1, v, std::dec)) fail = true; }
            else if (k == "maxNumReads") { if (!from_string<int>(maxNumReads,0, 200000, v, std::dec)) fail = true; }
            else if (k == "showCallReads") { if (!from_string<int>(showCallReads,0, 1, v, std::dec)) fail = true; }
            else if (k == "minNumHaplotypeOverlaps") { if (!from_string<int>(minNumHaplotypeOverlaps,0, 100, v, std::dec)) fail = true; }
            else if (k == "excludeSamplesFile") { if (!v.empty()) { excludeSamplesFile = v; } else { fail = true; } }
            else if (k == "maxNumCandidatesPerWindow") { if (!from_string<int>(maxNumCandidatesPerWindow,0, 64, v, std::dec)) fail = true; }
            else if (k == "maxNumReadsWindow") { if (!from_string<int>(maxNumReadsWindow,0, 100000, v, std::dec)) fail = true; }
            else if (k == "singleSampleHetThreshold") { if (!from_string<double>(singleSampleHetThreshold,0, 1000, v, std::dec)) fail = true; }
            else if (k == "singleSampleHomThreshold") { if (!from_string<double>(singleSampleHomThreshold,0, 1000, v, std::dec)) fail = true; }
            else if (k == "EMtol") { if (!from_string<double>(EMtol,0, 1, v, std::dec)) fail = true; }
            else if (k == "EMmaxiter") { if (!from_string<int>(EMmaxiter,0, 10000, v, std::dec)) fail = true; }
            else if (k == "doEM") { if (!from_string<int>(doEM,0, 1, v, std::dec)) fail = true; }
            else throw std::string("Unrecognized parameter: "+k);

            if (fail) throw std::string("Cannot determine value for parameter " + k + " from "+v);
      }
        else throw std::string("Incorrectly formatted parameter string: two consecutive colons");
    }

    checkAndInit();
}






/*


 VCF File


 */


std::vector<std::string> VCFFile::VCFEntry::getFilters() const
{
    std::vector<std::string> vFilters;
    Maps::const_iterator it = filters.begin();
    for (;it!=filters.end();it++) vFilters.push_back(it->first);
    return vFilters;
}

std::string VCFFile::VCFEntry::getInfoValue(const std::string& key) const
{
    Maps::const_iterator it = info.find(key);
    if (it == info.end()) return "-"; else return it->second;
}

void VCFFile::VCFEntry::parseInfoString()
{
    std::vector<std::string> infos = SplitString(this->infoString,';');

    for (size_t x=0;x<infos.size();x++)
    {
        std::vector<std::string> kv = SplitString(infos[x],'=');
        if (kv[0].empty()) throw std::string("VCFEntry: cannot parse info tag "+kv[0]);
        if (kv.size()==1) info[kv[0]]="TAG";
        else if (kv.size()==2) info[kv[0]]=kv[1];
        else throw std::string("VCFEntry: cannot parse info string");
    }

}

void VCFFile::VCFEntry::setFilters(const std::string & str)
{
     std::vector<std::string> f = SplitString(str,';');
     for (size_t x=0;x<f.size();x++) filters[f[x]] = "PRESENT";
}


template<class T> bool VCFFile::VCFEntry::fromInfoTag(T & val, const std::string & tag, std::ios_base& (*f)(std::ios_base&))
{
    std::string tagv = getInfoValue(tag);
    if (tagv=="-") return false;
    if (!from_string(val, tagv, f)) std::cerr << "Candidate VCF File: could not parse INFO tag " << tag << std::endl;
    return true;
}

void VCFFile::VCFEntry::write(std::ostream & out) const
{
    out << chrom << " " << ref << " " << alt << " " << infoString << std::endl;

}

VCFFile::VCFFile()
{
    m_isOpen = false;
    m_mode = "null";
}

VCFFile::VCFFile(const std::string& fileName, const std::string & mode)
{
    m_isOpen = false;
    m_fileName = fileName;
    assert(mode == "w" || mode == "r");

    // open the file

    if (mode == "w")
    {
        m_outputFileHandle.open(fileName.c_str());
    if (m_outputFileHandle.is_open()) m_isOpen = true;

    }
    else if (mode == "r")
    {
        m_inputFileHandle.open(fileName.c_str());
    if (m_inputFileHandle.is_open()) m_isOpen = true;

    }

    if (!m_isOpen)
    {
        std::string msg = "VCFFile::VCFFile: Cannot open file: '";
        msg.append(fileName);
        msg.append("'");
        if(fileName.empty()) msg.append("EMPTY FILENAME");
        throw msg;
    }
    m_mode = mode;
}

void VCFFile::setSamples(const std::vector<std::string> & samples)
{
    assert(!samples.empty());
    m_samples = samples;
}

void VCFFile::outputHeader(const std::string & refFile, const std::string & paramString)
{
    assert(m_isOpen && m_mode == "w");
    m_outputFileHandle << "##fileformat=VCFv4.0" << std::endl;
    m_outputFileHandle << "##source=Dindel2.0" << std::endl;
    m_outputFileHandle << "##reference=" << refFile << std::endl;

    //m_outputFileHandle << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"<< std::endl;
    //m_outputFileHandle << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=AC,Number=1,Type=Float,Description=\"Allele count\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads preferentially aligning to variant haplotype on forward strand\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=NF,Number=1,Type=Integer,Description=\"Number of reads preferentially aligning to variant haplotype on reverse strand\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=HSR,Number=1,Type=Integer,Description=\"Number of haplotype-specific (informative) reads\">"<< std::endl;
    m_outputFileHandle << "##INFO=<ID=HPLen,Number=1,Type=Integer,Description=\"Homopolymer length\">"<< std::endl;
    //m_outputFileHandle << "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"<< std::endl;
    //m_outputFileHandle << "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"<< std::endl;
    //m_outputFileHandle << "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"<< std::endl;
    m_outputFileHandle << "##FILTER=<ID=Candidate,Description=\"Variant is candidate for testing\">"<< std::endl;
    m_outputFileHandle << "##FILTER=<ID=q10,Description=\"Quality below 10\">"<< std::endl;
    m_outputFileHandle << "##FILTER=<ID=NoCall,Description=\"No call made\">"<< std::endl;
    m_outputFileHandle << "##FILTER=<ID=LowQuality,Description=\"Low quality variant\">"<< std::endl;
    //m_outputFileHandle << "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">"<< std::endl;
    m_outputFileHandle << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"<< std::endl;
    m_outputFileHandle << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"<< std::endl;
    //m_outputFileHandle << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"<< std::endl;
    m_outputFileHandle << "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"<< std::endl;
    m_outputFileHandle << "##DindelParameters=\"" << paramString << "\"" << std::endl;

    
    m_outputFileHandle << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    if (m_samples.size()>0)
    {
        m_outputFileHandle << "\tFORMAT";
        for (size_t x=0;x<m_samples.size();x++) m_outputFileHandle << "\t" << m_samples[x];
    }

    m_outputFileHandle << std::endl;
}

VCFFile::VCFEntry VCFFile::getNextEntry()
{
    assert(m_mode=="r");


    if (m_inputFileHandle.eof()) return VCFFile::VCFEntry();

    std::string line;

    bool foundEntry = false;
    while (!foundEntry)
    {
    if (m_inputFileHandle.eof()) return VCFFile::VCFEntry();
    std::getline(m_inputFileHandle,line);
    if (line.size()<1) return VCFFile::VCFEntry();
    if (line[0]!='#') foundEntry=true; else line.clear();
    }
    
    if (DINDEL_DEBUG) std::cerr << "here" << std::endl;


    std::istringstream is(line);

    VCFEntry e;

    int idx = 0;
    while (!is.eof())
    {
        std::string field;
        is >> field;
        switch(idx)
        {
            case 0: e.chrom=field; break;
            case 1:
                if (!from_string<int>(e.pos, field, std::dec)) throw std::string("VCF: cannot convert position in entry "+line);
                break;
            case 2: e.id = field; break;
            case 3: e.ref = field; break;
            case 4: e.alt = field; break;
            case 5:
                float qual;
                if (!from_string<float>(qual, field, std::dec)) throw std::string("VCF: cannot convert quality in entry "+line);
                e.qual = int(qual);
                break;
            case 6: e.setFilters(field); break;
            case 7: e.infoString = field; e.parseInfoString(); break;
        }
        idx++;
    }
    e.m_isEmpty = false;
    
    if (idx<7) throw std::string("VCF: line did not have all required fields " + line);

    if (DINDEL_DEBUG)
    {
    std::cerr << "getNextEntry(): ";
    e.write(std::cerr);
    }

    return e;
}

VCFFile::~VCFFile()
{
   if (m_mode == "w")
    m_outputFileHandle.close();
   else if (m_mode == "r")
    m_inputFileHandle.close();
}





DindelVariant addInsertion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, const std::string & seq)
{
    // insertion is added BEFORE refPos (see extractIndelsFromBamLeftJustified why)

    // determine left and right boundaries

    const std::string & refSeq = refHap.getSequence();
    int refseqlen = refSeq.length();

    int relPos = refPos-refHap.getRefStart();
    if (relPos<=0) throw std::string("relpos error");

    int len = (int) seq.length();

    int left=relPos-REPOSITION_INDEL_WINDOW;
    //while (left>=0 && refSeq.substr(left,len) == seq) left-=len;
    if (left<0) left=0;

    int right=relPos+REPOSITION_INDEL_WINDOW;
    if (right>refseqlen) right=refseqlen;
    //while (right+len<=refseqlen && refSeq.substr(right,len)==seq) right+=len;

    std::string eir = refSeq.substr(left, right-left+1);
    std::string alt = eir;
    alt.insert(relPos-left, seq);

    if (DINDEL_DEBUG)
    {

        std::cerr << "\n ***** INS \n ";
        std::cerr << "refPos: " << refPos << " relPos: " << relPos << " left: " << left << " right: " << right << " len: " << len << " seq: " << seq << "\n";
        std::cerr << refSeq.substr(relPos-25,25) << "|" << refSeq.substr(relPos,25) << "\n";
        std::cerr << "eir: " << eir << "\n";
        std::cerr << "alt: " << alt << "\n";
    }


    // shift to the left as far as possible

    int pos=relPos-left;
    std::string newSeq=seq;
    for (int p=relPos-left-1;p>=0;p--)
    {
        std::string _alt = eir;
        std::string nseq =alt.substr(p,len);
       _alt.insert(p, nseq);
       if (_alt == alt) {
           pos=p;
           newSeq=nseq;
       }

    }
    int newPos = pos+left-1; // -1 is for conversion to VCF4 convention

    std::string ref; ref += refSeq[newPos];
    alt=ref+newSeq;

    if (DINDEL_DEBUG)
    {
    std::cerr << " ref: " << ref << " alt: " << alt << " pos: " << newPos << std::endl;

    }
    return DindelVariant(chrom, ref, alt, newPos+refHap.getRefStart());
    
}

DindelVariant addDeletion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, int len)
{
    const std::string & refSeq = refHap.getSequence();
    int refseqlen = refSeq.length();
    int relPos = refPos-refHap.getRefStart();
    assert(relPos>0);

    assert(relPos+len<=int(refSeq.length()));
    std::string seq = refSeq.substr(relPos, len);

    int left=relPos-REPOSITION_INDEL_WINDOW;
    //while (left>=0 && refSeq.substr(left,len) == seq) left-=len;
    if (left<0) left=0;

    int right=relPos+REPOSITION_INDEL_WINDOW;
    if (right>refseqlen) right=refseqlen;
 
    std::string eir = refSeq.substr(left, right-left+1);
    std::string alt = eir;
    alt.erase(relPos-left, len);

    if (DINDEL_DEBUG)
    {

        std::cerr << "\n ***** DEL \n ";
        std::cerr << "refPos: " << refPos << " relPos: " << relPos << " left: " << left << " right: " << right << " len: " << len << "\n";
        std::cerr << refSeq.substr(relPos-25,25) << "|" << refSeq.substr(relPos,25) << "\n";
        std::cerr << "eir: " << eir << "\n";
        std::cerr << "alt: " << alt << "\n";
    }

    // shift to the left as far as possible

    int pos=relPos-left;
    std::string newSeq=eir.substr(pos, len);
    for (int p=relPos-left-1;p>=0;p--)
    {
        std::string _alt = eir;
        std::string nseq =eir.substr(p,len);
       _alt.erase(p, nseq.length());
       if (_alt == alt) {
           pos=p;
           newSeq=nseq;
       }

    }
    int newPos = pos+left-1;

    std::string ref; ref += refSeq[newPos]; 
    alt=ref;
    ref+=newSeq;
    assert(newSeq[0]==refSeq[newPos+1]);
    if (DINDEL_DEBUG)
    {
        std::cerr << "ref: " << ref << " alt: " << alt << " newPos: " << newPos << std::endl;
        std::cerr <<"\n";
    }
    return DindelVariant(chrom, ref, alt, newPos+refHap.getRefStart());

}




