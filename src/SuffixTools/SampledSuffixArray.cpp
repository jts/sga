//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SampledSuffixArray - Data structure holding
// a subset of suffix array positions. From the
// sampled positions, the other entries of the
// suffix array can be calculated.
//
#include "SampledSuffixArray.h"
#include "SAReader.h"

//
SampledSuffixArray::SampledSuffixArray()
{

}

SampledSuffixArray::SampledSuffixArray(const std::string& filename)
{
    // Read the sampled suffix array from a file
    
}

// 
SAElem SampledSuffixArray::calcSA(int64_t idx, const BWT* pBWT) const
{
    size_t offset = 0;
    SAElem elem;

    while(1)
    {
        if(idx % m_sampleRate == 0 && !m_saSamples[idx / m_sampleRate].isEmpty())
        {
            // A valid sample is stored for this idx
            elem = m_saSamples[idx / m_sampleRate];
            break;
        }

        // A sample does not exist for this position, perform a backtracking step
        char b = pBWT->getChar(idx);
        idx = pBWT->getPC(b) + pBWT->getOcc(b, idx - 1);

        if(b == '$')
        {
            // idx (before the update) corresponds to the start of a read.
            // We can directly look up the saElem for idx from the lexicographic index
            assert(idx < (int64_t)m_saLexoIndex.size());
            elem = m_saLexoIndex[idx];
            break;
        }
        else
        {
            // A backtracking step is performed, increment offset
            offset += 1;
        }
    }

    elem.setPos(elem.getPos() + offset);
    return elem;
}

// 
void SampledSuffixArray::build(std::string bwtFilename, std::string saiFilename, int sampleRate)
{
    m_sampleRate = sampleRate;

    // Read in the sampled suffix array
    SAReader saReader(saiFilename);

    size_t numStrings = 0;
    size_t numSAIElems = 0;
    saReader.readHeader(numStrings, numSAIElems);

    m_saLexoIndex.reserve(numStrings);
    saReader.readElems(m_saLexoIndex);

    // Load BWT
    BWT* pBWT = new BWT(bwtFilename);

    // Set the size of the sampled vector
    size_t numElems = (pBWT->getBWLen() / m_sampleRate) + 1;

    m_saSamples.resize(numElems);

    // For each sample, backtrack through the BWT until a valid sample is found
    for(size_t i = 0; i < numElems; ++i)
    {
        size_t idx = i * m_sampleRate;
        assert(idx < pBWT->getBWLen());
        SAElem elem = calcSA(idx, pBWT);
        m_saSamples[i] = elem;
    }

    delete pBWT;
}

// Validate the sampled suffix array values are correct
void SampledSuffixArray::validate(const std::string filename)
{
    ReadTable* pRT = new ReadTable(filename);
    SuffixArray* pSA = new SuffixArray(pRT, 1);
    BWT* pBWT = new BWT(pSA, pRT);
    
    std::cout << "Validating sampled suffix array entries\n";

    for(size_t i = 0; i < pSA->getSize(); ++i)
    {
        SAElem calc = calcSA(i, pBWT);
        SAElem real = pSA->get(i);
        if(calc.getID() != real.getID() || calc.getPos() != real.getPos())
        {
            std::cout << "Error SA elements do not match for " << i << "\n";
            std::cout << "Calc: " << calc << "\n";
            std::cout << "Real: " << real << "\n";
            exit(1);
        }
    }
    
    std::cout << "All calculate SA values are correct\n";

    delete pRT;
    delete pSA;
    delete pBWT;
}


