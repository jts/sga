//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapBlock - Data structures holding
// the result of the alignment of a sequence read
// to a BWT
// 
#ifndef OVERLAPBLOCK_H
#define OVERLAPBLOCK_H

#include "BWTInterval.h"
#include "SearchHistory.h"
#include "BitChar.h"
#include "BWT.h"
#include "SearchHistory.h"
#include "GraphCommon.h"

// Flags indicating how a given read was aligned to the FM-index
// Used for internal bookkeeping
struct AlignFlags
{
    public:
        AlignFlags() {}
        AlignFlags(bool qr, bool tr, bool qc)
        {
            setQueryRev(qr);
            setTargetRev(tr);
            setQueryComp(qc);
        }

        bool isQueryRev() const { return data.test(QUERYREV_BIT); }
        bool isTargetRev() const { return data.test(TARGETREV_BIT);    }
        bool isQueryComp() const { return data.test(QUERYCOMP_BIT); }

        // Returns true if the relationship between the query and target is reverse complement
        bool isReverseComplement() const { return isTargetRev() != isQueryRev(); }

        friend std::ostream& operator<<(std::ostream& out, const AlignFlags& af)
        {
            out << af.data;
            return out;
        }

        friend std::istream& operator>>(std::istream& in, AlignFlags& af)
        {
            in >> af.data;
            return in;
        }

        void write(std::ostream& out)
        {
            data.write(out);
        }

        void read(std::istream& in)
        {
            data.read(in);
        }

    private:

        void setQueryRev(bool b) { data.set(QUERYREV_BIT, b); }
        void setTargetRev(bool b) { data.set(TARGETREV_BIT, b); }
        void setQueryComp(bool b) { data.set(QUERYCOMP_BIT, b); }

        static const uint8_t QUERYREV_BIT = 0;
        static const uint8_t TARGETREV_BIT = 1;
        static const uint8_t QUERYCOMP_BIT = 2;
        BitChar data;
};


// A BWTInterval pair and associated alignment data
struct OverlapBlock
{
    OverlapBlock() {}
    
    OverlapBlock(BWTIntervalPair r, int ol, 
                 int nd, const AlignFlags& af)  : ranges(r), 
                                                  overlapLen(ol), 
                                                  numDiff(nd),
                                                  flags(af), isEliminated(false) {}

    OverlapBlock(BWTIntervalPair r, int ol, 
                 int nd, const AlignFlags& af, 
                 const SearchHistoryVector& backHist);

    // Returns the string that corresponds to this overlap block.
    // This is constructed by transforming the original string using the back
    // history while correcting for the reverse-complement searches
    std::string getOverlapString(const std::string& original);

    // Get the full string that is indicated by this overlap
    // The overlap cannot be a containment and it must have a 
    // forward extension history. The returned string is the string
    // of the actual read that forms the overlap with original. In other words
    // it might not be the same strand as the original read.
    std::string getFullString(const std::string& original);

    // Return a pointer to the BWT that should be used to extend the block
    // this is the opposite BWT that was used in the backwards search
    const BWT* getExtensionBWT(const BWT* pBWT, const BWT* pRevBWT) const;

    // Return the spectrum of extensions given by the interval in ranges
    // The counts are given in the canonical frame, which means that
    // if the query string was reversed, we flip the counts
    AlphaCount getCanonicalExtCount(const BWT* pBWT, const BWT* pRevBWT) const;

    // Return the index of the interval corresponding to the frame of 
    // reference for the original read
    int getCanonicalIntervalIndex() const;

    // Return the direction of the edge that this overlap block describes
    EdgeDir getEdgeDir() const;

    // Construct an overlap record from this overlap block
    Overlap toOverlap(const std::string queryID, const std::string targetID, int queryLen, int targetLen) const;

    // Comparison operator, compare by lower coordinate of 0 
    friend bool operator<(const OverlapBlock& a, const OverlapBlock& b)
    {
        return a.ranges.interval[0].lower < b.ranges.interval[0].lower;    
    }

    static bool sortSizeDescending(const OverlapBlock& ob1, const OverlapBlock& ob2)
    {
        return ob1.overlapLen > ob2.overlapLen;
    }

    static bool sortIntervalLeft(const OverlapBlock& ob1, const OverlapBlock& ob2)
    {
        return ob1.ranges.interval[0].lower < ob2.ranges.interval[0].lower;
    }

    // I/O
    friend std::ostream& operator<<(std::ostream& out, const OverlapBlock& obl)
    {
        out << obl.ranges << " " << obl.overlapLen << " " << obl.numDiff << " " << obl.flags;
        return out;
    }

    friend std::istream& operator>>(std::istream& in, OverlapBlock& obl)
    {
        in >> obl.ranges >> obl.overlapLen >> obl.numDiff >> obl.flags;
        return in;
    }

    // Data
    BWTIntervalPair ranges;
    int overlapLen;
    int numDiff;
    AlignFlags flags;
    bool isEliminated;

    // The sequence of divergent bases during the backward and forward search steps
    SearchHistoryVector backHistory;
    SearchHistoryVector forwardHistory;
};

// Collections
typedef std::list<OverlapBlock> OverlapBlockList;
typedef OverlapBlockList::iterator OBLIter;

// Global Functions

// Ensure all the overlap blocks in the list are distinct
void removeSubMaximalBlocks(OverlapBlockList* pList);

// Given the overlapping blocks A and B, construct a list of blocks where the index ranges do not intersect
OverlapBlockList resolveOverlap(const OverlapBlock& A, const OverlapBlock& B);

// Partition the overlap block list into two lists, 
// one for the containment overlaps and one for the proper overlaps
void partitionBlockList(int readLen, OverlapBlockList* pCompleteList, 
                        OverlapBlockList* pOverlapList, 
                        OverlapBlockList* pContainList);

#endif
