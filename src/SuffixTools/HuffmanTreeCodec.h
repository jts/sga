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
    size_t code;
    size_t bits;
};

// Huffman Tree Implementation
template<typename T>
class HuffmanTreeCodec
{
    public:
       
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
        struct DecodePair
        {
            T symbol;
            uint8_t bits;
        };

        //
        struct HuffmanNodePriority
        {
            bool operator()(HuffmanNode* pLHS, HuffmanNode* pRHS) const
            {
                return (pLHS->frequency > pRHS->frequency);
            }
        };

        struct HuffmanTriplet
        {
            static bool compareCodes(const HuffmanTriplet& lhs, const HuffmanTriplet& rhs)
            {
                return lhs.code < rhs.code;
            }

            T symbol;
            size_t code;
            size_t bits;
        };
        
        // typedefs
        typedef std::map<T, int> CountMap;
        typedef std::map<T, EncodePair> EncodeTable;
        typedef std::vector<DecodePair> DecodeTable;
        typedef std::vector<HuffmanNode*> HuffmanNodePtrVector;
        typedef std::vector<HuffmanTriplet> TripletVector;
        typedef std::priority_queue<HuffmanNode*,
                                    std::vector<HuffmanNode*>, 
                                    HuffmanNodePriority > HuffmanQueue;

    

        //
        HuffmanTreeCodec()
        {

        }

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
            buildEncodeTable(huffQueue.top(), 0, 0, "");
            buildDecodeTable();

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

        DecodePair decode(size_t code) const
        {
            assert(code < m_decoder.size());
            return m_decoder[code];
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

        void buildEncodeTable(HuffmanNode* pNode, size_t currCode, size_t codeBits, std::string codeStr)
        {
            if(pNode->isLeaf())
            {
                //std::cout << "HC " << pNode->symbol << " f: " << pNode->frequency << " " << currCode << " " << codeBits << " " << codeStr << "\n";
                if(m_minSymbolBits == 0 || codeBits < m_minSymbolBits)
                    m_minSymbolBits = codeBits;

                if(codeBits > m_maxSymbolBits)
                    m_maxSymbolBits = codeBits;

                EncodePair ep = {currCode, codeBits};
                m_encoder.insert(std::make_pair(pNode->symbol, ep));
            }
            else
            {
                if(pNode->pLeftChild)
                    buildEncodeTable(pNode->pLeftChild, (currCode << 1), codeBits + 1, codeStr + "0");
                if(pNode->pRightChild)
                    buildEncodeTable(pNode->pRightChild, (currCode << 1) + 1, codeBits + 1, codeStr + "1");
            }
        }
        
        // Construct a decoder table
        void buildDecodeTable()
        {
            assert(!m_encoder.empty());

            // Construct Huffman Triplets representing the symbol encodings
            TripletVector tripletVector;
            for(typename EncodeTable::iterator iter = m_encoder.begin(); iter != m_encoder.end(); ++iter)
            {
                HuffmanTriplet triplet = {iter->first, iter->second.code, iter->second.bits};
                tripletVector.push_back(triplet);
            }

            // Sort triplets by code
            std::sort(tripletVector.begin(), tripletVector.end(), HuffmanTriplet::compareCodes);

            //size_t max = tripletVector.back().code;
            size_t maxBits = m_maxSymbolBits;
            size_t maxCode = (1 << maxBits) - 1;
            m_decoder.reserve(maxCode);

            for(size_t i = 0; i <= maxCode; ++i)
            {
                size_t idx = findPrefixIdx(i, maxBits, tripletVector);
                DecodePair dp;
                if(idx != tripletVector.size())
                {
                    dp.symbol = tripletVector[idx].symbol;
                    dp.bits = tripletVector[idx].bits;
                }
                else
                {
                    dp.bits = 0; // signal the code isnt used
                }
                //std::cout << "DT: " << int2Binary(i, 8) << " " << (int)dp.symbol << " " << (int)dp.bits << "\n";
                m_decoder.push_back(dp);
            }
        }

        // Search a triplet vector for the index of a code that is a prefix of code y
        // m is the number of bits in code y
        // If no search code exists, return vec.size()
        size_t findPrefixIdx(size_t y, size_t m, TripletVector& vec)
        {
            size_t i = 0;
            while(i <= vec.size())
            {
                if(isPrefix(vec[i].code, y, vec[i].bits, m))
                    return i;
                ++i;
            }
            return i;
        }
        
        // Return true if x is a prefix of y where x is a n bit string and y is an m bit string
        bool isPrefix(size_t x, size_t y, size_t n, size_t m)
        {
            // shift x until it is the same number of bits as m
            x <<= (m - n);
            // Construct a mask where the first n bits out of m are set
            size_t mask = 0;

            // Set the low n bits
            for(size_t i = 0; i < n; ++i)
            {
                mask <<= 1;
                mask += 1;
            }

            // Shift into place
            mask <<= (m - n);
            return (x & mask) == (y & mask);
        }

        size_t m_minSymbolBits;
        size_t m_maxSymbolBits;
        EncodeTable m_encoder;    
        DecodeTable m_decoder;
};

#endif
