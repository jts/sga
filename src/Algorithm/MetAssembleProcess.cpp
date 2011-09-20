///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MetAssembleProcess - Find contigs in an
// abstractly represented de Bruijn graph
// of a metagenome
//
#include "MetAssembleProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "MetagenomeBuilder.h"

//
//
//
MetAssemble::MetAssemble(const MetAssembleParameters& params) : m_parameters(params)
{
}

//
MetAssemble::~MetAssemble()
{
}

MetAssembleResult MetAssemble::process(const SequenceWorkItem& item)
{
    MetAssembleResult result;
    SeqRecord currRead = item.read;
    std::string w = item.read.seq.toString();
    if(w.size() <= m_parameters.kmer)
    {
        return result;
    }

    // Perform a backwards search using the read sequence
    // Check which k-mers have already been visited using the
    // shared bitvector. If any bit in the range [l,u] is set
    // for a suffix of the read, then we do not visit those kmers
    // later
    int len = w.size();
    int num_kmers = len - m_parameters.kmer + 1;
    std::vector<bool> visitedKmers(false, num_kmers);

    int j = len - 1;
    char curr = w[j];
    BWTInterval interval;
    BWTAlgorithms::initInterval(interval, curr, m_parameters.pBWT);
    --j;

    for(;j >= 0; --j)
    {
        curr = w[j];
        BWTAlgorithms::updateInterval(interval, curr, m_parameters.pBWT);
        assert(interval.isValid());

        // At this point interval represents the suffix [j,len)
        // Check if the starting point of this interval is set
        if(j < num_kmers)
            visitedKmers[j] = m_parameters.pBitVector->test(interval.lower);
    }
    
    // Process the kmers that have not been previously visited
    for(j = 0; j < num_kmers; ++j)
    {
        if(visitedKmers[j])
            continue; // skip

        std::string kmer = w.substr(j, m_parameters.kmer);
        
        // Get the interval for this kmer
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, kmer);

        // Check if this interval has been marked by a previous iteration of the loop
        assert(interval.isValid());
        if(m_parameters.pBitVector->test(interval.lower))
            continue;

        BWTInterval rc_interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, reverseComplement(kmer));

        size_t count = interval.size();
        if(rc_interval.isValid())
            count += rc_interval.size();

        if(count >= m_parameters.kmerThreshold)
        {
            // Process the kmer
            std::string contig = processKmer(kmer, count);

            // We must determine if this contig has been assembled by another thread.
            // Break the contig into lexicographically ordered set of kmers. The lowest kmer is chosen
            // to represent the contig. If this kmer has been marked as visited, we discard the contig
            // otherwise we mark all kmers in the contig and output the contig.
            StringVector kmers = getLexicographicKmers(contig);

            // Get the lowest kmer in the set
            assert(!kmers.empty());
            std::string lowest = kmers.front();
            BWTInterval lowInterval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, lowest);
            BWTInterval lowRCInterval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, reverseComplement(lowest));
            bool marked = false;

            // If the kmer exists in the read set (the interval is valid), we attempt to mark it
            // Otherwise, we attempt to mark its reverse complement. If either call succeeds
            // we output the contig. This block of code gives the synchronization between threads.
            // If multiple threads attempt to assemble the same contig, marked will be true for only
            // one of the threads.
            if(lowInterval.isValid())
                marked = m_parameters.pBitVector->updateCAS(lowInterval.lower, false, true);
            else
                marked = m_parameters.pBitVector->updateCAS(lowRCInterval.lower, false, true);
            markSequenceKmers(contig);
            if(marked)
            {
                if(contig.size() >= m_parameters.minLength)
                    result.contigs.push_back(contig);
            }
        }

        // Update the bit vector for the source kmer
        for(int64_t i = interval.lower; i <= interval.upper; ++i)
            m_parameters.pBitVector->updateCAS(i, false, true);

        for(int64_t i = rc_interval.lower; i <= rc_interval.upper; ++i)
            m_parameters.pBitVector->updateCAS(i, false, true);
    }
        
    return result;
}

//
std::string MetAssemble::processKmer(const std::string& str, int count)
{
    MetagenomeBuilder builder;
    builder.setSource(str, count);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT, m_parameters.pBWTCache, m_parameters.pRevBWTCache);
    builder.run();

    StringVector contigs;
    builder.getContigs(contigs);
    std::string out = contigs.empty() ? "" : contigs.front();
    //std::cout << "Constructed " << contigs.size() << " contigs from kmer (size: " << out.size() << ")\n";
    return out;
}

// Update the bit vector with the kmers that were assembled into str
void MetAssemble::markSequenceKmers(const std::string& str)
{
    assert(str.size() >= m_parameters.kmer);
    size_t n = str.size() - m_parameters.kmer + 1;

    for(size_t i = 0; i < n; ++i)
    {
        std::string kseq = str.substr(i, m_parameters.kmer);
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_parameters.pBitVector->updateCAS(j, false, true);
        }

        // Mark the reverse complement k-mers too
        std::string rc_kseq = reverseComplement(kseq);
        interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, rc_kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_parameters.pBitVector->updateCAS(j, false, true);
        }
    }
}

// Break a contig into a lexicographically ordered set of kmers
StringVector MetAssemble::getLexicographicKmers(const std::string& contig) const
{
    StringVector out;
    for(size_t i = 0; i < contig.size() - m_parameters.kmer + 1; ++i)
    {
        std::string ks = contig.substr(i, m_parameters.kmer);
        std::string rs = reverseComplement(ks);
        out.push_back(ks < rs ? ks : rs);
    }

    std::sort(out.begin(), out.end());
    return out;
}

//
// MetAssembleAggregateResult
//
MetAssembleAggregateResults::MetAssembleAggregateResults(const std::string& filename) : m_numContigs(0), m_basesWritten(0)
{
    //
    m_pWriter = createWriter(filename);
}

//
MetAssembleAggregateResults::~MetAssembleAggregateResults()
{
    delete m_pWriter;
}

void MetAssembleAggregateResults::process(const SequenceWorkItem& /*item*/, const MetAssembleResult& result)
{
    // Write to the contigs file
    for(size_t i = 0; i < result.contigs.size(); ++i)
    {
        std::stringstream idMaker;
        idMaker << "contig-" << m_numContigs++;

        SeqItem item = { idMaker.str(), result.contigs[i] };
        m_basesWritten += result.contigs[i].length();
        item.write(*m_pWriter);
    }
}
