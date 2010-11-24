//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// HuffmanTreeCodec -- templated implementation of a huffman tree
//
#ifndef HUFFMANTREECODEC_H
#define HUFFMANTREECODEC_H

#include <inttypes.h>
#include <map>
#include <vector>
#include <string>
#include <queue>

struct EncodePair
{
    uint8_t code;
    uint8_t bits;
};

struct DecodePair
{
    char base;
    uint8_t bits;
};

// Huffman Tree Implementation
template<typename T>
class HuffmanTreeCodec
{
    // Node used in the tree
    struct HuffmanNode
    {
        HuffmanNode() : pParent(NULL), pLeftChild(NULL), pRightChild(NULL), isInternal(false), frequency(0.0f) {}
        bool isLeaf() const { return pLeftChild == NULL && pRightChild == NULL; }

        // Data
        HuffmanNode* pParent;
        HuffmanNode* pLeftChild;
        HuffmanNode* pRightChild;

        T symbol;
        bool isInternal;
        double frequency;
    };
    
    //
    struct HuffmanNodePriority
    {
        bool operator()(HuffmanNode* pLHS, HuffmanNode* pRHS) const
        {
            return (pLHS->frequency > pRHS->frequency);
        }
    };
    
    // typedefs
    typedef std::map<T, int> CountMap;
    typedef std::map<T, EncodePair> EncodeTable;
    typedef std::vector<HuffmanNode*> HuffmanNodePtrVector;
    typedef std::priority_queue<HuffmanNode*,
                                std::vector<HuffmanNode*>, 
                                HuffmanNodePriority > HuffmanQueue;

    public:
    

        // Construct the tree using a map of symbol counts
        HuffmanTreeCodec(const CountMap& input)
        {
            //
            size_t sum = 0;
            for(typename CountMap::const_iterator iter = input.begin(); iter != input.end(); ++iter)
            {
                sum += iter->second;
            }

            // Create initial leaves and place in priority queue
            HuffmanQueue huffQueue;
            HuffmanNodePtrVector nodePtrVector;
            for(typename CountMap::const_iterator iter = input.begin(); iter != input.end(); ++iter)
            {
                HuffmanNode* pNode = new HuffmanNode;
                pNode->symbol = iter->first;
                pNode->frequency =  (double)iter->second / sum;
                huffQueue.push(pNode);
                nodePtrVector.push_back(pNode);
            }

            while(huffQueue.size() > 1)
            {
                HuffmanNode* pLowest = huffQueue.top();
                huffQueue.pop();
                HuffmanNode* pSecond = huffQueue.top();
                huffQueue.pop();

                //std::cout << "Lowest: " << pLowest->symbol << ", " << pLowest->frequency << "\n";
                //std::cout << "Second: " << pSecond->symbol << ", " << pSecond->frequency << "\n";

                // Merge these nodes
                double new_freq = pLowest->frequency + pSecond->frequency;
                HuffmanNode* pMerged = new HuffmanNode;
                pMerged->frequency = new_freq;
                pMerged->pLeftChild = pSecond;
                pMerged->pRightChild = pLowest;
                huffQueue.push(pMerged);
                nodePtrVector.push_back(pMerged);
            }
            
            // Traverse the tree building the codes
            m_minSymbolBits = 0;
            m_maxSymbolBits = 0;
            buildEncodeTable(huffQueue.top(), "");

            // Delete the tree
            for(size_t i = 0; i < nodePtrVector.size(); ++i)
            {
                delete nodePtrVector[i];
                nodePtrVector[i] = 0;
            }
        }

        EncodePair encode(T sym) const
        {
            typename EncodeTable::const_iterator iter = m_encoder.find(sym);
            assert(iter != m_encoder.end());
            return iter->second;
        }

        size_t getMinBits() const { return m_minSymbolBits; }
        size_t getMaxBits() const { return m_maxSymbolBits; }

        // Returns the greatest symbol in the set of encoding values
        // whose value is no greater than val
        T getGreatestLowerBound(T val) const
        {
            typename EncodeTable::const_iterator iter = m_encoder.upper_bound(val);
            assert(iter == m_encoder.end() || iter->first > val);
            assert(iter != m_encoder.begin());
            iter--;
            T lower = iter->first;
            assert(lower <= val);
            return lower;
        }

    private:

        void buildEncodeTable(HuffmanNode* pNode, std::string code)
        {
            if(pNode->isLeaf())
            {
                //std::cout << "Huff " << pNode->symbol << " f: " << pNode->frequency << " " << code.size() << " " << code << "\n";
                if(m_minSymbolBits == 0 || code.size() < m_minSymbolBits)
                    m_minSymbolBits = code.size();

                if(code.size() > m_maxSymbolBits)
                    m_maxSymbolBits = code.size();

                EncodePair ep = {0, code.size()};
                m_encoder.insert(std::make_pair(pNode->symbol, ep));
            }
            else
            {
                if(pNode->pLeftChild)
                    buildEncodeTable(pNode->pLeftChild, code + "0");
                if(pNode->pRightChild)
                    buildEncodeTable(pNode->pRightChild, code + "1");
            }
        }

    private:
        size_t m_minSymbolBits;
        size_t m_maxSymbolBits;
        EncodeTable m_encoder;    

};

#endif
