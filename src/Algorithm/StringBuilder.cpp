///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringBuilder - Iteratively construct strings
// that represent sequences in an assembly graph.
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "StringBuilder.h"
#include "BWTAlgorithms.h"
#include "LRAlignment.h"
#include "StdAlnTools.h"

//
// StringThreaderNode
//
StringThreaderNode::StringThreaderNode(const std::string& l, 
                                       int tae, int qae, int as) : label(l),
                                                                   target_alignment_end(tae),
                                                                   query_alignment_end(qae),
                                                                   alignment_score(as)
{


}

// Delete the children of the node
StringThreaderNode::~StringThreaderNode()
{
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;
}

//
void StringThreaderNode::printFullAlignment(const std::string* pQuery) const
{
    StdAlnTools::printGlobalAlignment(label, *pQuery);
}

//
// StringTheader
//
StringThreader::StringThreader(const std::string& seed, 
                               const std::string* pQuery, 
                               int kmer, 
                               const BWT* pBWT, 
                               const BWT* pRevBWT) : m_pBWT(pBWT), m_pRevBWT(pRevBWT), 
                                                     m_kmer(kmer), m_pQuery(pQuery)
{
    // Create the root node containing the seed string
    WARN_ONCE("TODO: FIX INITIAL ALIGNMENT COORDINATES");
    m_pRootNode = new StringThreaderNode(seed, seed.size(), seed.size(), 0);
}

//
StringThreader::~StringThreader()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

// Run the threading algorithm
void StringThreader::run()
{
    m_pRootNode->printFullAlignment(m_pQuery);
}

//
// String Builder. Deprecated
//
StringBuilder::StringBuilder(const std::string& seed, int kmer, const BWT* pBWT, const BWT* pRevBWT) : m_pBWT(pBWT), m_pRevBWT(pRevBWT), m_kmer(kmer)
{
    m_strings.push_back(seed);
}

// Perform one round of extension for all strings in the collection
void StringBuilder::extendOnce()
{
    StringVector newStrings;
    for(size_t i = 0; i < m_strings.size(); ++i)
    {
        // Get the spectrum of extensions for this string
        std::string& curr = m_strings[i];

        // do not try to extend deleted strings
        if(curr.empty())
            continue;
        
        // Get the last k-1 bases of this string
        std::string s = curr.substr(curr.size() - m_kmer + 1, m_kmer - 1);
        BWTIntervalPair intervalPair = BWTAlgorithms::findIntervalPair(m_pBWT, m_pRevBWT, s);
        
        std::cout << s << ": "<< intervalPair << "\n";
        assert(intervalPair.interval[0].isValid());

        AlphaCount64 extCount = BWTAlgorithms::getExtCount(intervalPair.interval[1], m_pRevBWT);
        int num_extension = 0;
        char first_extension = '\0';

        for(size_t j = 0; j < DNA_ALPHABET::size; ++j)
        {
            char b = DNA_ALPHABET::getBase(j);
            if(extCount.get(b) > 0)
            {
                // Perform extension or branch
                if(num_extension == 0)
                {
                    // Defer the extension of curr until any branches have been formed
                    first_extension = b;
                }
                else
                {
                    newStrings.push_back(curr + b);
                }
                num_extension += 1;
            }
        }
        
        if(num_extension > 0)
            curr.append(1, first_extension);
    }

    m_strings.insert(m_strings.end(), newStrings.begin(), newStrings.end());
}

// Cull low-scoring strings from the collection
void StringBuilder::cull(const std::string& query, size_t n)
{
    if(m_strings.size() <= n)
        return;

    //    
    // Calculate scores between all the strings in the vector and the query
    //

    // Set up alignment parameters
    LRAlignment::LRParams params;
    AlnParam par;
    int matrix[25];
    StdAlnTools::setAlnParam(par, matrix, params.alnParams);

    // Make a packed version of the query
    IntVector scores;

    uint8_t* pQueryT = StdAlnTools::createPacked(query);
    for(size_t i = 0; i < m_strings.size(); ++i)
    {
        // Skip deleted strings
        const std::string& target = m_strings[i];
        if(target.empty())
            continue;

        uint8_t* pTargetT = StdAlnTools::createPacked(target);
        path_t path;
        int score = aln_extend_core(pTargetT, target.size(), pQueryT, query.size(), &par, &path, 0, 1, NULL);
        
        //int tal = path.i;
        //int qal = path.j;
        //printf("Score: %d tl: %zu tal: %d qal: %d\n", score, target.size(), tal, qal);
        scores.push_back(score);
        delete [] pTargetT;
    }
    delete [] pQueryT;

    // Choose the score cutoff
    IntVector sortedScores = scores;
    std::sort(sortedScores.begin(), sortedScores.end(), std::greater<int>());
    int cutoff = sortedScores[n];

    // Keep only strings with score at least cutoff. This may be larger than cutoff size.
    StringVector tempStrings;
    for(size_t i = 0; i < m_strings.size(); ++i)
    {
        if(scores[i] >= cutoff)
            tempStrings.push_back(m_strings[i]);
    }
    m_strings.swap(tempStrings);

    //std::cout << m_strings.size() << " strings after cull (cutoff: " << cutoff << "\n";
}

// Print the set of strings found
void StringBuilder::print() const
{
    std::cout << "SB has " << m_strings.size() << "\n";
    for(size_t i = 0; i < m_strings.size(); ++i)
    {
        std::cout << i << ": " << m_strings[i] << "\n";
    }
}

// Print the pairwise alignment of a query string against the set
void StringBuilder::print(const std::string& query) const
{
    //
    // Set up alignment parameters
    //
    LRAlignment::LRParams params;
	size_t max_target = StdAlnTools::calculateMaxTargetLength(query.size(), params.alnParams);
    int max_path_length = max_target + query.size();
    path_t* path = (path_t*)calloc(max_path_length, sizeof(path_t));    

    AlnParam par;
    int matrix[25];
    StdAlnTools::setAlnParam(par, matrix, params.alnParams);

    // Make a packed version of the query
    uint8_t* pQueryT = StdAlnTools::createPacked(query);

    MAlignDataVector mAlignVec;
    for(size_t i = 0; i < m_strings.size(); ++i)
    {
        std::string target = m_strings[i];
        if(target.empty()) // deleted string
            continue;
        
        //
        // Perform the alignment
        //
        uint8_t* pTargetT = StdAlnTools::createPacked(target);
        int path_len = 0;
		int score = aln_global_core(pTargetT, target.size(), pQueryT, query.size(), &par, path, &path_len);
        assert(path_len <= max_path_length);

        int cigarLen = 0;
		uint32_t* pCigar = aln_path2cigar32(path, path_len, &cigarLen);
        
        //
        MAlignData maData;
        maData.str = target;
        maData.position = 0;
        maData.targetID = i;
        maData.targetAlignLength = score; //target.size();
        maData.targetPosition = -1;

        // Add initial padding to cigar
        maData.expandedCigar.append(maData.position, 'S');
        
        // Add alignment symbols
        for (int j = 0; j != cigarLen; ++j)
            maData.expandedCigar.append(pCigar[j]>>4, "MID"[pCigar[j]&0xf]);

        mAlignVec.push_back(maData);
        delete [] pTargetT;
    }

    MultiAlignment ma(query, mAlignVec);
    ma.print();
    delete [] pQueryT;
    free(path);
}
