//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#include <algorithm>
#include <iostream>
#include "MultiOverlap.h"
#include "Alphabet.h"

MultiOverlap::MultiOverlap(const std::string& rootID, const std::string& rootSeq) : m_rootID(rootID), m_rootSeq(rootSeq)
{
    
}

//
void MultiOverlap::add(const std::string& seq, const Overlap& ovr)
{
    MOData mod(seq, ovr);

    // Swap root read into first position if necessary
    if(mod.ovr.id[0] != m_rootID)
        mod.ovr.swap();
    assert(mod.ovr.id[0] == m_rootID);

    // RC the sequence if it is different orientation than the root
    if(mod.ovr.match.isRC())
    {
        mod.seq = reverseComplement(mod.seq);
        mod.ovr.match.canonize();
    }

    // Initialize the offset value, the amount that a coordinate 
    // for the non-root sequence must be shifted so that
    // the sequences are aligned
    mod.offset = mod.ovr.match.inverseTranslate(0);
    m_overlaps.push_back(mod);
}

//
void MultiOverlap::add(const MOData& mod)
{
    assert(mod.ovr.id[0] == m_rootID);
    m_overlaps.push_back(mod);
}

//
Overlap MultiOverlap::getOverlap(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].ovr;
}

//
std::string MultiOverlap::simpleConsensus() const
{
    std::string out;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount ac = p.getAlphaCount();
        char maxBase;
        BaseCount maxCount;
        ac.getMax(maxBase, maxCount);
        BaseCount rootCount = ac.get(m_rootSeq[i]);

        if(rootCount == maxCount)
            out.push_back(m_rootSeq[i]);
        else
            out.push_back(maxBase);
    }
    return out;
}

// Determine if this multi-overlap has a conflicted position,
// a position in the multioverlap where the second-most 
// prevalent base has a frequency greater than cutoff
bool MultiOverlap::isConflicted(size_t cutoff) const
{
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount ac = p.getAlphaCount();

        char order[5];
        ac.getSorted(order, 5);
    
        // If the second-most prevalent base is greater than the cutoff
        // we consider the MO to be conflicted
        if(ac.get(order[1]) > cutoff)
            return true;
    }
    return false;
}

// Returns true if the overlap at idx has the include flag set
int MultiOverlap::getPartition(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].partitionID;
}

// Returns true if the overlap at idx has the include flag set
void MultiOverlap::setPartition(size_t idx, int p)
{
    assert(idx < m_overlaps.size());
    m_overlaps[idx].partitionID = p;
}

// Return the total number of bases in the multioverlap
size_t MultiOverlap::getNumBases() const
{
    size_t count = 0;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        count += p.getDepth();
    }
    return count;
}

//
std::string MultiOverlap::consensusConflict(double p_error, int conflictCutoff)
{
    std::vector<AlphaCount> acVec;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount ac = p.getAlphaCount();
        acVec.push_back(ac);
    }

    for(size_t j = 0; j < m_overlaps.size(); ++j)
    {
        std::string cstr;
        std::string sc;
        std::string nc;

        int score = 0;
        int sum = 0;
        int wrong = 0;
        for(size_t i = 0; i < acVec.size(); ++i)
        {
            char rootBase = m_rootSeq[i];
            char sorted[ALPHABET_SIZE];
            acVec[i].getSorted(sorted, ALPHABET_SIZE);
            int second = acVec[i].get(sorted[1]);
            bool isConflict = second > conflictCutoff;
            char b = getMODBase(m_overlaps[j], i);

            if(isConflict)
            {
                if(b != '\0' && (rootBase == sorted[0] || rootBase == sorted[1]))
                {
                    if(b == rootBase)
                        ++score;
                    else
                        ++wrong;
                    ++sum;
                }

                if(b == sorted[0] || b == sorted[1])
                {
                    nc.push_back(b);
                }
                else
                {
                    nc.push_back('N');
                }
            }
            
            if(b != '\0')
                cstr.push_back(b);    
            else
                cstr.push_back('N');
        }

        double frac;
        if(sum == 0)
            frac = 1.0f;
        else
            frac = (double)score/(double)sum;

        m_overlaps[j].score = frac;
        if(sum > 0 && wrong > 0)
            m_overlaps[j].partitionID = 1;
        else
            m_overlaps[j].partitionID = 0;

        //printf("CFT\t*%s\t%d,%lf\t%s\n", cstr.c_str(), (int)j, frac, nc.c_str());
    }

    /*
    std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);
    for(size_t i = 0; i < m_overlaps.size() && i < 15; ++i)
    {
        m_overlaps[i].partitionID = 0;
    }
    */
    std::string consensus = calculateConsensusFromPartition(p_error);
    for(size_t i = 0; i < acVec.size(); ++i)
    {
        char sorted[ALPHABET_SIZE];
        acVec[i].getSorted(sorted, ALPHABET_SIZE);
        int second = acVec[i].get(sorted[1]);
        bool isConflict = second > conflictCutoff;
        if(!isConflict)
        {
            consensus[i] = sorted[0];
        }
    }
    
    return consensus;
}

void MultiOverlap::partitionBest(double p_error, size_t n)
{
    // Compute the likelihood of the alignment between the root sequence
    // and every member of the multioverlap
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        double likelihood = 0.0f;
        double total_bases = 0;
        double mismatches = 0;
        double overlap_len = 0;
        for(size_t j = 0; j < m_rootSeq.size(); ++j)
        {
            Pileup pileup = getSingletonPileup(j, i);
            if(pileup.getDepth() == 2)
            {
                ++overlap_len;
                if(pileup.getBase(0) != pileup.getBase(1))
                {
                    ++mismatches;
                }
            }
            DNADouble ap = pileup.calculateLikelihoodNoQuality(p_error);
            likelihood += ap.marginalize(0.25f);
            total_bases += pileup.getDepth();
        }

        m_overlaps[i].score = likelihood / total_bases;
        m_overlaps[i].partitionID = 1;
        //std::cout << "ERS\t" << m_overlaps[i].score << "\t" << likelihood << "\n";
    }

    std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);

    for(size_t i = 0; i < m_overlaps.size() && i < n; ++i)
        m_overlaps[i].partitionID = 0;
}

// Construct the left and right seqtries representing the multioverlap
void MultiOverlap::makeSeqTries(double p_error, SeqTrie& leftTrie, SeqTrie& rightTrie)
{
    double lp = log(p_error);
    if(m_overlaps.empty())
        return;

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        MOData data = m_overlaps[i];

        if(data.offset == 0)
            leftTrie.insert(data.seq, lp);
        else
            rightTrie.insert(reverse(data.seq), lp);
    }    
}

//
size_t MultiOverlap::countPartition(int id) const
{
    size_t count = 0;
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        if(m_overlaps[i].partitionID == id)
            ++count;
    }
    return count;
}

std::string MultiOverlap::calculateConsensusFromPartition(double p_error)
{
    std::string out;
    out.reserve(m_rootSeq.size());

    // require the best base call to be above this above to correct it
    double epsilon = 0.01;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p0;
        Pileup p1;
        getPartitionedPileup(i, p0, p1);
        DNADouble ap = p0.calculateLikelihoodNoQuality(p_error);
        char best_c;
        char curr_c = m_rootSeq[i];
        double max;
        ap.getMax(best_c, max);

        // Require the called value to be substantially better than the
        // current base
        if(best_c == curr_c || (best_c != curr_c && max - ap.get(curr_c) < epsilon))
        {
            out.push_back(curr_c);
        }
        else
        {
            out.push_back(best_c);
        }
    }
    return out;
}

// Compute the likelihood of the multioverlap
double MultiOverlap::calculateLikelihood() const
{
    WARN_ONCE("Compute likelihood using fixed error rate");

    double likelihood = 0.0f;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup pileup = getPileup(i);
        DNADouble ap = pileup.calculateLikelihoodNoQuality(0.01);
        likelihood += ap.marginalize(0.25f);
    }
    return likelihood;
}

// Calculate the likelihood of the multioverlap given the groups defined by IntVector
// There are only two possible groups, 0 and 1. The root_sequence (index 0 in the vector) is assumed to belong
// to group 1
double MultiOverlap::calculateGroupedLikelihood() const
{
    double likelihood = 0.0f;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup g0;
        Pileup g1;
        getPartitionedPileup(i, g0, g1);
        if(g0.getDepth() > 0)
            likelihood += g0.calculateLikelihoodNoQuality(0.01).marginalize(0.25f);
        if(g1.getDepth() > 0)
            likelihood += g1.calculateLikelihoodNoQuality(0.01).marginalize(0.25f);
    }
    return likelihood;
}


// Calculate the probability of the 4 bases given the multi-overlap
// for a single position
DNADouble MultiOverlap::calcAlphaProb(size_t idx) const
{
    Pileup pileup = getPileup(idx);
    return pileup.calculateSimpleAlphaProb();
}

// Get the number of times each base appears at position
// idx in the multi-align
AlphaCount MultiOverlap::calcAlphaCount(size_t idx) const
{
    Pileup pileup = getPileup(idx);
    return pileup.getAlphaCount();
}


// Calculate the probability of this multioverlap
void MultiOverlap::calcProb() const
{
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup pileup = getPileup(i);
        if(pileup.getDepth() > 1)
        {
            char cnsBase = pileup.calculateSimpleConsensus();
            DNADouble ap = pileup.calculateSimpleAlphaProb();
            char refBase = pileup.getBase(0);
            
            if(refBase != cnsBase)
            {
                int ref_count = pileup.getCount(pileup.getBase(0));
                int cons_count = pileup.getCount(cnsBase);
                int depth = pileup.getDepth();
                double rp = ap.get(refBase);
                double cp = ap.get(cnsBase);
                printf("CNS\t%d\t%d\t%d\t%lf\t%lf\n", ref_count, cons_count, depth, rp, cp);
            }
        }
    }

    /*
    int numMismatches = 0;
    int numAlignedBases = 0;
    double errorRate = 0.01;

    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup pileup = getPileup(i);
        if(pileup.getDepth() > 1)
        {
            char consensus = pileup.calculateSimpleConsensus();
            
            // Calculate the number of bases in the pileup that do not match the consensus
            for(size_t j = 0; j < pileup.getDepth(); ++j)
            {
                if(pileup.getBase(j) != consensus)
                    ++numMismatches;
                ++numAlignedBases;
            }
        }
    }

    double actualRate = double(numMismatches) / double(numAlignedBases);
    double expectedMismatches = errorRate * double(numAlignedBases);
    printf("MM\t%d\t%lf\t%d\t%lf\t%lf\n", numMismatches, expectedMismatches, numAlignedBases, actualRate, errorRate);
    */
}

// Return the base in mod that matches the base at
// idx in the root seq. If mod does not overlap 
// the root at this position, returns '\0'
char MultiOverlap::getMODBase(const MOData& mod, int idx) const
{
    int trans_idx = idx - mod.offset;
    if(trans_idx >= 0 && size_t(trans_idx) < mod.seq.size())
        return mod.seq[trans_idx];
    else
        return '\0';
}

// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base
Pileup MultiOverlap::getPileup(int idx) const
{
    Pileup p;
    p.add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            p.add(curr.seq[trans_idx]);
        }
    }
    return p;
}

// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base only including the first numElems items
Pileup MultiOverlap::getPileup(int idx, int numElems) const
{
    assert((size_t)numElems < m_overlaps.size());
    Pileup p;
    p.add(m_rootSeq[idx]);

    for(int i = 0; i < numElems; ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && (size_t)trans_idx < curr.seq.size())
        {
            p.add(curr.seq[trans_idx]);
        }
    }
    return p;
}


// Get the pileup between the root sequence and one of the sequences in the MO
Pileup MultiOverlap::getSingletonPileup(int base_idx, int ovr_idx) const
{
    assert(ovr_idx < (int)m_overlaps.size());

    Pileup p;
    p.add(m_rootSeq[base_idx]);

    const MOData& curr = m_overlaps[ovr_idx];
    // translate idx into the frame of the current sequence
    int trans_idx = base_idx - curr.offset;
    if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
    {
        p.add(curr.seq[trans_idx]);
    }
    return p;
}

// Fill in the pileup groups g0 and g1
void MultiOverlap::getPartitionedPileup(int idx, Pileup& g0, Pileup& g1) const
{
    g0.add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            if(curr.partitionID == 0)
                g0.add(curr.seq[trans_idx]);
            else
                g1.add(curr.seq[trans_idx]);
        }
    }
}


// Fill in the pileup groups g0 and g1
PileupVector MultiOverlap::getPartitionedPileup(int idx, int num_parts) const
{
    assert(false && "untested");
    PileupVector pv(num_parts);
    pv.push_back(Pileup());
    
    // the root always goes in partition 0
    pv.back().add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            assert(curr.partitionID < num_parts);
            pv[curr.partitionID].add(curr.seq[trans_idx]);
        }
    }
    return pv;
}

// Print the MultiOverlap groups specified by the IntVec to stdout
void MultiOverlap::printGroups()
{
    for(int i = 0; i <= 1; ++i)
    {
        MultiOverlap mo_group(m_rootID, m_rootSeq);
        for(size_t j = 0; j < m_overlaps.size(); ++j)
        {
            const MOData& curr = m_overlaps[j];
            if(curr.partitionID == i)
                mo_group.add(curr);
        }
        std::cout << "MO GROUP " << i << "\n";
        mo_group.print();
    }
}

// Print the MultiOverlap to stdout
void MultiOverlap::print(int default_padding, int max_overhang)
{
    std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::sortOffset);
    std::cout << "\nDrawing overlaps for read " << m_rootID << "\n";
    int root_len = int(m_rootSeq.size());
    
    // Print the root row at the bottom
    printRow(default_padding, max_overhang, root_len, 0, root_len, 0, 0.0f, m_rootSeq, m_rootID);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        int overlap_len = curr.ovr.match.getMaxOverlapLength();
        int nd = curr.ovr.match.countDifferences(m_rootSeq, curr.seq);
        double score = static_cast<double>(nd) / overlap_len;
        printRow(default_padding, max_overhang, root_len, 
                 curr.offset, overlap_len, nd, 
                 score, curr.seq, curr.ovr.id[1].c_str());
    }    
}

// Print a single row of a multi-overlap to stdout
void MultiOverlap::printRow(int default_padding, int max_overhang, 
                            int root_len, int offset, int overlap_len, int pid, 
                            double score, const std::string& seq, const std::string& id)
{
    int c_len = seq.length();

    // This string runs from c_offset to c_offset + len
    // Clip the string at -max_overhang to root_len + max_overhang
    int left_clip = std::max(offset, -max_overhang);
    int right_clip = std::min(offset + c_len, root_len + max_overhang);
    
    // translate the clipping coordinates to the string coords
    int t_left_clip = left_clip - offset;
    int t_right_clip = right_clip - offset;
    
    // Calculate the length of the left padding
    int padding = default_padding + left_clip;
    std::string leader = (t_left_clip > 0) ? "..." : "";
    std::string trailer = (t_right_clip < c_len) ? "..." : ""; 
    std::string clipped = seq.substr(t_left_clip, t_right_clip - t_left_clip);
    padding -= leader.size();

    assert(padding >= 0);
    std::string padding_str(padding, ' ');
    std::string outstr = padding_str + leader + clipped + trailer;
    printf("%s\t%d\t%d\t%lf\tID:%s\n", outstr.c_str(), overlap_len, pid, score, id.c_str());
}

void MultiOverlap::printMasked()
{
    printf("ROOT\t%s\t%s\n", m_rootSeq.c_str(), m_rootID.c_str());

    for(size_t j = 0; j < m_overlaps.size(); ++j)
    {
        std::string out;
        for(size_t i = 0; i < m_rootSeq.length(); ++i)
        {
            char b = getMODBase(m_overlaps[j], i);
            if(b == '\0')
                out.push_back('.');
            else if(b == m_rootSeq[i])
                out.push_back('-');
            else
                out.push_back(b);
        }
        printf("OVRL\t%s\t%s\n", out.c_str(), m_overlaps[j].ovr.id[1].c_str());
    }
}

// Print the MultiOverlap horizontally, in a pileup format
void MultiOverlap::printPileup()
{
    std::cout << "\nDrawing overlap pileup for read " << m_rootID << "\n";
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        printf("%zu\t%s\n", i, p.toStr().c_str());
    }
}

//
bool MultiOverlap::MOData::sortOffset(const MOData& a, const MOData& b)
{
    return a.offset < b.offset;
}

// Sort by the ID of the non-root sequence, which is 
// guarenteed to bethe second id in the overlap structure
bool MultiOverlap::MOData::sortID(const MOData& a, const MOData& b)
{
    assert(a.ovr.id[0] == b.ovr.id[0]);
    return a.ovr.id[1] < b.ovr.id[1];
}

