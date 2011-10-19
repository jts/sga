//-----------------------------------------------
// Copyright 201- Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
//        and Kees Albers (caa@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DindelRealignWindow - Infer haplotypes using the FM-index
//

#ifndef DINDELREALIGNWINDOW_H
#define DINDELREALIGNWINDOW_H
#include <fstream>
#include "Util.h"
#include "HashMap.h"
#include "MultiAlignment.h"
#include <iomanip>
#include <list>
#include <set>
#include <map>


// Enums
enum DindelWindowCandHapAlgorithm {LINEAR} ;


// Constants
const size_t DINDEL_HASH_SIZE=8;
const int DINDEL_HMM_BANDWIDTH=6;

// Event types
const int DELETION_NOVEL_SEQUENCE = -1;
const int SNP = -2;
const int MULTINUCLEOTIDE_RUN = -3;
const int INSERTION = -4;
const int LEFTOVERHANG = -5;
const int RIGHTOVERHANG = -6;
const int ADDVARIANT_DEBUG = -9;

// Typedefs
typedef std::string SampleName;

// Forward declares
class DindelRealignWindow;

// Functions
long long int combinations(const int n, const int k);

// Complement a single base
inline char _complement(char base)
{
    switch(base)
    {
        case 'A':
            return 'T';
        case'C':
            return 'G';
        case 'G':
            return 'C';
        case'T':
            return 'A';
        default:
            return 'N';
    }
}

// Add the log-scaled values l1 and l2 using a transform to avoid
// precision errors
inline double addLogs(const double l1, const double l2)
{
    if (l1>l2) {
            double diff=l2-l1;
            return l1+log(1.0+exp(diff));
    } else {
            double diff=l1-l2;
            return l2+log(1.0+exp(diff));
    }
}

// Convert the type t into a string. Returns true if successful
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

// Convert the string s into a value of type T. Clamp
// the out value at min or max if the input string is out of the range
template <class T>
bool from_string(T& t, T min, T max,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  bool failed = !(iss >> f >> t).fail();

  if (!failed)
  {
      if (t<min) t=min;
      else if (t>max) t=max;
  }
  return failed;
}

// Split a string into tokens based on the given separator
// JS TODO: Util.h has a similar function, can deduplicate
std::vector<std::string> SplitString(const std::string & str, char sep);

// Parse a colon-separated string into a chromosome/start/end triple
void parseRegionString(const std::string & region, std::string & chrom, int & start, int & end);

//
// VariantPriors - Prior probability values for various variant types
//
class VariantPriors
{
    public:
        VariantPriors();
        double getDefaultProbVariant(const std::string & type) const;
        double m_probSNP, m_probINDEL, m_probMNP;

};

//
// DindelRead - A sequence read used as input to the dindel realignment
//
class DindelRead
{
    public:

        // Constructor
        DindelRead(const SeqItem & seqItem, const SampleName & sampleName, double mappingQual, int fixedBaseQual, bool isForward);

        // Functions

        const std::string  getSequence() const { return m_seqItem.seq.toString(); }
        const std::vector<double> getBaseQuals() const;
        const std::string getID() const { return m_seqItem.id; }
        double getMappingQual() const { return m_mappingQual; } //return double(bam->core.qual); }
        double getLogProbNotMapping() const { return m_mappingQual*0.23026; }
        bool isForward() const { return m_isForward; } //FIXME
        bool mateIsForward() const { assert(false); return false; } //FIXME
        
        //BAM bool BAMCigarHasIndel() const;
        //BAM int getBamStartPos() const { return bam->core.pos; }
        //BAM int getBamStartPosAdjusted() const;
        //BAM int getBAMEnd() const { return bam->core.n_cigar? bam_calend(&bam->core, bam1_cigar(bam)) : bam->core.pos + 1; }
        //BAM int getBAMEndAdjusted() const;
        
        bool isUnmapped() const {  return false; } //FIXME
        bool mateIsUnmapped() const { assert(false); return false; } //FIXME
        char getBase(int b) const { return m_seqItem.seq.get(size_t(b)); }
        int getQual(int b) const { b=1; return m_fixedBaseQual; }
        int length() const { return int(m_seqItem.seq.length()); }
        
        //BAM const bam1_t *getBam() const { return bam; }

        void setMappingQuality(double mappingQual) { m_mappingQual=mappingQual; }
        void setRCRead(bool status) { m_rcRead=status; if (status) setupHash(); } // need to recompute hash
        bool getRCRead() const { return m_rcRead; }
        void getLogProbCorrectError(std::vector<double> & lpCorrect, std::vector<double> & lpError) const;
        const std::vector<unsigned int> & getHashKeys();
        const std::string & getSampleName() const { return m_sampleName; }
        const std::string getLibraryName() const { return std::string("NAN"); }//FIXME

    private:
        
        // Functions
        void setupHash();

        // Data
        //BAM const bam1_t *bam;
        //BAM const DindelBAM *m_pDindelBAM;
        double m_mappingQual;
        int m_fixedBaseQual;
        bool m_rcRead, m_isForward, m_setupHash; // reverse-complement the read. Intended for unmapped reads. Determined using the status of the mate.
        SampleName m_sampleName;
        SeqItem m_seqItem;
        
        std::vector<unsigned int> m_hashKeys; // stores hashes for the read

};

//
// DindelVariant - A difference from the reference haplotype
//
class DindelVariant
{
    public:

        // Constructos
        DindelVariant();      
        DindelVariant(const std::string & chrom, const std::string & ref, const std::string & alt, int pos);

        // Functions
        const std::string & getChrom() const { return m_chrom; }
        const std::string & getRef() const { return m_ref; }
        const std::string & getAlt() const { return m_alt; }

        int getPos() const { return m_pos; }
        bool operator<(const DindelVariant & v) const { if (m_chrom<v.m_chrom) return true; else if (m_pos<v.m_pos) return true; else return m_id<v.m_id; }
        //bool operator<(const DindelVariant & v) const { return m_id < v.m_id; }

        // checks ref and pos against reference sequence
        //void check(Fasta *pFasta);
        // set min and max positions in haplotype where variant may be positioned without affecting the sequence
        void setHaplotypeUnique(int leftUniquePos, int rightUniquePos) { m_leftUniquePos = leftUniquePos; m_rightUniquePos = rightUniquePos; };
        int getHaplotypeLeftUnique() const;
        int getHaplotypeRightUnique() const;
        int getLengthDifference() const { return m_dlen; };
        //int getIndex() const;
        //void setIndex(int idx);
        std::string getID() const { return m_id; }
        static bool variantFromWindowString(const std::string & varstring, DindelVariant & variant);
        const std::string & getType() const { return m_type; }
        void setPriorProb(double prob) { m_priorProb = prob; }
        double getPriorProb() const { if (m_priorProb==-1.0) { std::cerr << "ERRVAR: "; write(std::cerr); } assert(m_priorProb!=-1.0); return m_priorProb; }
        void write(std::ostream & out) const;

    private:

        // Functions
        void setup();
        void checkSequence(const std::string & seq); // checks whether sequence has no funny characters

        // Data
        std::string m_chrom;
        std::string m_ref;
        std::string m_alt;

        std::string m_id;
        std::string m_type;

        double m_priorProb; // prior probability

        int m_pos;

        // how far can indel be shifted to left and right?
        int m_leftUniquePos, m_rightUniquePos;

        // does the variant change the length of the reference sequence?
        // m_dlen>0 is insertion
        int m_dlen;

        int m_idx; // index of variant in vector of DindelWindow


};

//
// DindelSequenceHash
//
class DindelSequenceHash
{
    public:
        
        // Constructors
        DindelSequenceHash() {};
        DindelSequenceHash(const std::string & sequence);

        // Functions
        const std::set<int> & lookup(const std::string & seq, int pos) const;
        //int lookup(unsigned int key) const;
        const std::list<int> & lookupList(unsigned int key) const;
        
        static inline unsigned int convert(const std::string & seq, int pos)
        {
            unsigned int v=0;
            for (int x=pos, y=0;x<pos+int(DINDEL_HASH_SIZE);x++,y++) {
                v |= (map_char(seq[x]) << (2*y) );
            }
            return v;
        }
    
        static inline unsigned int pushBack(const unsigned int & key, const char & c)
        {
            return (key >> 2) | (map_char(c) << (2*(DINDEL_HASH_SIZE-1)));
        }

        static inline int map_char(const char & c)
        {
            switch (c)
            {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T' : return 3;
            }
            return 0;
        }
 
        static void getKeys(std::vector<unsigned int> & keys, const std::string & sequence);
        void print() const;

    protected:
        typedef HashMap<unsigned int, std::list<int> > Hash;
        typedef HashMap<unsigned int, int > UniqueHash;

        // Functions
        void makeHash(const std::string & sequence);

        // Data
        Hash m_hash;
        // UniqueHash m_uniqueHash;
        std::set<int> emptySet;
        std::list<int> emptyList;
};

//
// DindelHaplotype
//
class DindelHaplotype
{
    public:

        // Constructors
        // a DindelHaplotype must be constructed starting from the reference haplotype.
        // Differences with the reference can be added by calling addVariant
        DindelHaplotype(const std::string & refSeq, int refSeqStart, bool isReference);
        DindelHaplotype(const std::string & refName, const std::string & refSeq, int refSeqStart, const MultiAlignment & ma, size_t varRow, size_t refRow);
        DindelHaplotype(const DindelHaplotype & haplotype, int copyOptions);

        // Functions
        bool addVariant(const DindelVariant & var);
        const std::vector<DindelVariant> getVariants() const { return m_variants; }
        const std::string & getSequence() const { return m_seq; }
        bool isReference() const { return m_isReference; } 
        void write(std::ostream & out) const;
        int getHomopolymerLength(int b) const { return m_hplen[b]; }
        int getHomopolymerLengthRefPos(int refPos) const { int b=getHapBase(refPos); assert(m_refPos[b]==refPos); return m_hplen[b]; }
        int length() const { return (int) m_seq.length();}
        int getRefBase(int hapBase) const { return m_refPos[hapBase]; };
        int getHapBase(int refPos) const;
        int getRefStart() const { return m_hapRefStart; }
        bool hasVariant(const DindelVariant & variant) const ;
        int getClosestDistance(const DindelVariant & variant, int hapPos1, int hapPos2) const;
        const DindelSequenceHash & getHash() const { return m_sequenceHash; }

    private:
        
        // Functions

        // initializes haplotype from reference sequence
        // sets homopolymer length vector
        void initHaplotype();
        void determineHomopolymerLengths();
        void updateHaplotype();
        
        // Data
        std::string m_seq;
        std::vector<int> m_hplen; // gives the homopolymer run length for any base in the haplotype sequence going left and right.
        std::vector<int> m_refPos; // position of haplotype base on reference.

        // list of VCF4 style variants contained in haplotype
        std::vector<DindelVariant> m_variants;
        HashMap<std::string, std::pair<int, int> > m_variant_to_pos;

        // position of leftmost haplotype base on reference sequence
        int m_hapRefStart;

        bool m_isReference;
    
        // constants for m_refPos
        DindelSequenceHash m_sequenceHash;
};

//
// DindelWindow contains a list of variants, the candidate haplotypes generated from those variants, and the relationships
// between the candidate haplotypes and the variants.
// The candidate haplotypes are annotated with respect to the reference sequence.
// The candidate variants are specified as modifications to the reference sequence by any string.replace operation
//
class DindelWindow
{
    public:
        
        // Constructor
        DindelWindow(const std::vector<DindelVariant> & variants, const std::string & refHap, int refHapStart);

        // Create window from a set of haplotypes and a reference sequence.
        // Uses SGA MultiAlignment to annotate the variations in the haplotypes with respect to the reference sequence.
        DindelWindow(const std::vector<std::string> & haplotypeSequences, const std::string & refHap, int refHapStart, const std::string & refName);

        ~DindelWindow();

        // Functions
  
        const std::vector<DindelHaplotype> & getHaplotypes() const { return haplotypes; }
        const std::string getChrom() const { return m_chrom; }
        int getLeft() const { return m_leftRef; }
        int getRight() const { return m_rightRef; }
        void addVariant(const DindelVariant & variant, bool addToAll);
        void filterHaplotypes(const std::vector<bool> & filterHaplotype);

    private:
    
        // Functions
        void makeWindow(const std::vector<DindelVariant> & variants);
        DindelHaplotype makeAltHapFromReference(const DindelVariant & var );
        void generateCandidateHaplotypes(const std::vector<DindelVariant> & variants);
        void setupRefSeq(const std::string & refHap, int refHapStart);
        bool addHaplotype(const DindelHaplotype & haplotype);
        void generateCandidateHaplotypesFromMultiAlignment();

        // Data
        
        // candidate haplotypes
        std::vector<DindelHaplotype> haplotypes; 
      
        // leftmost and rightmost coordinates of window on reference sequence
        int m_leftRef, m_rightRef;
        
        // reference sequence for [m_leftRef, m_rightRef]
        std::string m_refSeq;
        std::string m_chrom;

        // links variant to haplotype
        // std::vector<std::vector<int> > m_variantToHaplotype;

        // keeps track of haplotypes.
        std::map<std::string, int> m_hashAltHaps;

        //Fasta *m_pFasta; // pointer to Fasta class for getting the reference sequence

        //parameters
        DindelWindowCandHapAlgorithm m_candHapAlgorithm;
        MultiAlignment *m_pHaplotype_ma;
};

//
// VCFFile - Reader/Writer for VCF files
//
class VCFFile
{
    public:
        class VCFEntry
        {
            friend class VCFFile;
            public:
                
                // Typedefs
                typedef std::map<std::string, std::string> Maps;
                
                // Constructors
                VCFEntry() { m_isEmpty = true; };

                // Functions
                bool isEmpty() { return m_isEmpty; };
                std::string getInfoValue(const std::string & key) const;
                std::vector<std::string> getFilters() const;
                void parseInfoString();
                
                template<class T> 
                bool fromInfoTag(T & val, 
                                 const std::string & tag, 
                                 std::ios_base& (*f)(std::ios_base&));

                void setFilters(const std::string & str);
                void write(std::ostream & out);

                // Data
                std::string chrom, ref, alt, id, infoString;
                int pos, qual;
                Maps  info, filters;

            private:

                // Data
                bool m_isEmpty;
        };

        // Constructors
        VCFFile();
        VCFFile(const std::string & fileName, const std::string & mode);
        ~VCFFile();

        // Functions
        std::ostream & getOutputStream() { assert(m_mode=="w"); return m_outputFileHandle; }
        VCFEntry getNextEntry();
        void setSamples(const std::vector<std::string> & samples);
        const std::vector<std::string> & getSamples() const { return m_samples; }
        void outputHeader(const std::string & refFile, const std::string & paramString);
        
    private:

        // Data
        std::ofstream m_outputFileHandle;
        std::ifstream m_inputFileHandle;
        std::string m_fileName;
        bool m_isOpen;
        std::string m_mode;
        std::vector<std::string> m_samples;

};

//
// CoalescentResult
//
struct CoalescentResult
{
    // Constructors
    CoalescentResult() { qual = -1.0; mapFreq = -1.0; singleSampleCall = false; }
    
    CoalescentResult(double _qual, double _mapFreq)
    {
        qual = _qual;
        mapFreq = _mapFreq;
        singleSampleCall = false;
        singleSampleQual = -1.0;
    }

    CoalescentResult(double _qual, double _mapFreq, bool _isSingleSampleCall, double _singleSampleQual)
    {
        qual = _qual;
        mapFreq = _mapFreq;
        singleSampleCall = _isSingleSampleCall;
        singleSampleQual = _singleSampleQual;
    }

    // Data
    double qual, singleSampleQual;
    double mapFreq; // MAP frequency of haplotype
    bool singleSampleCall; // call based on simply thresholding individual sample genotype/haplotype likelihoods.
};

//
// DindelRealignWindowResult
//
class DindelRealignWindowResult
{
    public:

        //
        class Inference
        {
            public:

                // Functions
                void outputAsVCF(const DindelVariant & var, 
                                 const DindelRealignWindowResult & result, 
                                 VCFFile & vcfFile) const;

                static double computeStrandBias(int numForward, int numReverse);
                
                // these functions will take care of scaling and initialization
                void addDistanceToHistogram(int distance);
                void addAlignLikToHistogram(double logLik); //
                void addMapQToHistogram(double mappingQuality);

                // Data

                // Phred-scaled posterior prob
                double qual;
                double freq;
                double strandBias;
                // std::set<int> readsForward, readsForwardZeroMismatch;
                // std::set<int> readsReverse, readsReverseZeroMismatch;
                int numReadsForward;
                int numReadsReverse;
                int numReadsForwardZeroMismatch;
                int numReadsReverseZeroMismatch;
                int numUnmapped;
                int numLibraries;
                int numReadNames;
                
                bool isSingleSampleCall;
                double singleSampleQual;

                std::vector<int> histDistance; // histogram for distances to variant
                std::vector<int> histAlignLik; // histogram for read-haplotype alignment likelihoods
                std::vector<int> histMapQ; // histogram for mapping qualities
                std::string infoStr; // special info string that can be set freely
        };

        class GenotypeCall
        {
            public:

                // Constructor
                GenotypeCall() { qual=0.0; count = 0; called = false;}

                // Data
                double qual;
                int count;
                bool called; // variant was present in called pair of haplotypes
                double gl[3];
        };

        // Typedefs
        typedef std::map<DindelVariant, Inference> VarToInference;
        typedef std::map<DindelVariant, GenotypeCall> VarToGenotypeCall;
        typedef HashMap<std::string, VarToGenotypeCall> SampleToGenotypes;

        // Constructors
        DindelRealignWindowResult() : m_pDindelRealignWindow(NULL){ };
        DindelRealignWindowResult(const DindelRealignWindow & dindelRealignWindow) : m_pDindelRealignWindow(&dindelRealignWindow){ };
    
        // Functions
        void outputVCF(VCFFile & vcfFile);

        // Data
        std::vector<DindelHaplotype> haplotypes;
        std::vector<double> haplotypeFrequencies;
        
        // integrates same variants across different haplotypes
        VarToInference variantInference; 
        SampleToGenotypes sampleToGenotypes;
        
        // Number of reads that were used to call variants 
        // (that map well to one haplotype but not another)
        int numHapSpecificReads; 

        // number of inferred haplotypes
        int numHaplotypes; 
        bool outputGenotypes;
        
    private:
        
        // Data
        const DindelRealignWindow *m_pDindelRealignWindow;
};

//
// CandidateSNP
//
class CandidateSNP
{
    public:
        // Constructor
        CandidateSNP(char refBase)
        {
            m_refBase = refBase;
            baseToReads = BaseToReads(4);
        }

        // Funtions
        char getBase(size_t idx) const { return "ACGT"[idx]; }
        size_t getBaseIdx(char base) const 
        { 
            switch(base) 
            { 
                case 'A' : return 0; 
                case 'C' : return 1; 
                case 'G' : return 2; 
                case 'T' : return 3; 
                case 'N' : throw std::string("N not allowed"); 
            }; 
            return 0; 
        }

        
        void addRead(char refBase, char altBase, size_t readIndex)
        {
            if (m_refBase!=refBase)
            {
                std::cerr << "ERROR: m_refBase: " << m_refBase << " refBase: " << refBase << std::endl;
            }
            assert(m_refBase==refBase);
            baseToReads[getBaseIdx(altBase)].insert(readIndex);
        }

        // Data
        typedef std::vector < std::set<size_t> > BaseToReads;
        BaseToReads baseToReads;
        char m_refBase;
};

//
// PosToCandidateSNP
//
class PosToCandidateSNP : public HashMap<int, CandidateSNP >
{
    public:
        void addSNP(int readIndex, int refPos, char refBase, char altBase);
};

//
// ReadHaplotypeAlignment
//
class ReadHaplotypeAlignment
{
    public:

        // Constructors
        ReadHaplotypeAlignment() 
        { 
            logLik = 1.0; 
            nmm =-1; 
            postProbLastReadBase=-1.0; 
            isUngapped=false; 
            hapPosLastReadBase=-1;
        }

        ReadHaplotypeAlignment(double _logLik, int _nmm) 
        { 
            logLik = _logLik; 
            nmm = _nmm; 
            postProbLastReadBase=-1.0; 
            isUngapped=false; 
            hapPosLastReadBase=-1;
        }
        
        ReadHaplotypeAlignment(double _logLik, int _nmm, bool _isUngapped) 
        { 
            logLik = _logLik; 
            nmm = _nmm; 
            postProbLastReadBase=-1.0; 
            isUngapped=_isUngapped; 
            hapPosLastReadBase=-1;
        }
       
        // Data
        double logLik;
        
        // posterior probability of last read base
        double postProbLastReadBase; 
        
        // number of mismatches between read and haplotype
        int nmm; 

        // position of last readbase on haplotype. -1 indicates off haplotype (shouldn't happen)
        int hapPosLastReadBase; 

        bool isUngapped;
};

//
// DindelRealignParameters
//
class DindelRealignParameters
{
    friend class DindelRealignAlgorithm;
    friend class DindelRealignWindow; 

    public:

        // Constructors
        DindelRealignParameters();
        DindelRealignParameters(const std::string & paramString);

    private:
        
        // Functions
        std::string getDefaultParameters() const;
        void initFromString(const std::string & paramString);
        const std::string & getParamString() const { return m_paramString; }
        void checkAndInit();

        // Data
        std::string m_paramString;
        std::string excludeSamplesFile;

        int windowReadBuffer;
        int minVariantSep, haplotypeWidth, minCandidateAlleleCount;
        int maxNumReads; // maximum number for fetch
        int maxNumReadsWindow; // maximum number of reads in window
        int maxNumCandidatesPerWindow;
        double theta_snp;
        double theta_indel;

        VariantPriors variantPriors;
        int maxMappingQuality;

        int addSNPMaxMismatches; // maximum number of mismatches allowed in ungapped alignment for adding candidate SNPs. Note this would include the SNP!
        double addSNPMinMappingQual; // minimum mapping quality of read for adding candidate SNPs
        int addSNPMinBaseQual; // minimum base quality
        size_t addSNPMinCount;
        int addSNPMaxSNPs;
        double minPostProbLastReadBaseForUngapped; // minimum posterior probability of the last base in order to attempt an ungapped alignment for SNP detection

        int genotyping; // genotyping mode? If 1, then each read will be aligned against each haplotype using the HMM

        int showCallReads; // show reads driving the call
        int minNumHaplotypeOverlaps; // in filterHaplotypes, how many reads should overlap the sequence surrounding the variant?

        // derived parameters
        double minLogLikAlignToRef, minLogLikAlignToAlt; // minimum logLik for aligning to a reference and an alternative haplotype
        double addSNPMinLogBaseQual;

        double singleSampleHetThreshold, singleSampleHomThreshold;

        int doEM;
        double EMtol;
        int EMmaxiter;
};

// Realigns reads in a given window.
// pOverlapper and dindelReads should correspond 1-1. It is assumed that pOverlapper gives back indices into the dindelReads array.
class DindelRealignWindow
{
    public:
        
        // Constructors
        DindelRealignWindow(const DindelWindow *pDindelWindow, 
                            std::vector<DindelRead> & dindelReads, 
                            const DindelRealignParameters & realignParameters);
        
        // Functions
        void run(const std::string & algorithm, VCFFile & vcfFile);

        const DindelWindow & getDindelWindow() const { return  m_dindelWindow; }
        const std::vector< std::vector< ReadHaplotypeAlignment > > & getHapReadAlignments() const 
        { 
            return hapReadAlignments; 
        }

    private:

        // Data
        DindelWindow m_dindelWindow;
        std::vector<DindelRead> *m_pDindelReads;
        DindelRealignParameters realignParameters;
        std::vector<int> m_readMapsToReference;

        // HAPLOTYPE ALIGNMENT BUSINESS
        PosToCandidateSNP m_posToCandidateSNP;

        // Functions
        void filterHaplotypes();
        void addSNPsToHaplotypes();
        void processHaplotypes(size_t firstHapIdx, size_t lastHapIdx);

        // convert seed positions into a likelihood given the candidate haplotype
        ReadHaplotypeAlignment computeReadHaplotypeAlignment(size_t readIndex, 
                                                             const DindelHaplotype& haplotype, 
                                                             int start, 
                                                             int end, 
                                                             const std::vector<double> & lpError, 
                                                             const std::vector<double> & lpCorrect, 
                                                             bool rcReadSeq);

        void HMMAlignReadAgainstHaplotypes(size_t readIndex, 
                                           size_t firstHap, 
                                           size_t lastHap, 
                                           const std::vector<double> & lpCorrect, 
                                           const std::vector<double> & lpError);

        // HAPLOTYPE FREQUENCY ESTIMATION BUSINESS

        // haplotype frequencies
        class HL
        {
            public:
                HL(double _ll, int _idx) { ll = _ll; idx = _idx;}

                bool operator<(const HL & hl) const 
                { 
                    if (int(ll*10.0)==int(hl.ll*10.0))
                    {
                        if (idx<hl.idx) 
                            return true; 
                        else 
                            return false;
                    } 
                    else
                    {
                        return int(ll*10.0)<int(hl.ll*10.0);
                    }
                }

                double ll;
                int idx;
        };

        // Functions
        void printReadAlignments(int readIdx, std::ostream & out, int offset, bool supportAlt);

        //
        void doEM(const std::vector< std::vector<double> > & hrLik, 
                  const std::vector<int> & calledHaplotypes, 
                  std::vector<double> & haplotypeFrequencies);

        //
        DindelRealignWindowResult estimateHaplotypeFrequencies(double minLogLikAlignToRef, 
                                                               double minLogLikAlignToAlt, 
                                                               bool capUsingMappingQuality, 
                                                               bool print);
        
        DindelRealignWindowResult estimateHaplotypeFrequenciesModelSelection(double minLogLikAlignToRef, 
                                                                             double minLogLikAlignToAlt, 
                                                                             bool capUsingMappingQuality, 
                                                                             bool print);

        double getHaplotypePrior(const DindelHaplotype & h1, const DindelHaplotype & h2) const;

        // this is a function that hides how the read information is actually stored in the various ReadTables
        const DindelRead & getRead(size_t readIndex) const { return (*m_pDindelReads)[readIndex]; }

        // haplotype-read likelihoods
        
        // first dim is haplotypes, second dim is reads
        std::vector< std::vector< ReadHaplotypeAlignment > > hapReadAlignments; 

        // update hapReadAlignments.
        void computeReadHaplotypeAlignmentsUsingHMM(size_t firstHap, size_t lastHap);
        void addDiploidGenotypes(DindelRealignWindowResult & result, bool useEstimatedHaplotypeFrequencies);

        typedef HashMap<int, CoalescentResult> HapIdxToCoalescentResult;

        class ReadSamples
        {
            public:
        
                //
                ReadSamples(const std::vector<DindelRead> & dindelReads);
                
                // Data
                std::vector<std::string> samples;
                std::vector<int> readIndices;
                std::vector<double> logNumStates;
            private:

                const std::vector<DindelRead> *m_pDindelReads;
        };

        ReadSamples m_readSamples;

        HapIdxToCoalescentResult computeHaplotypesCoalescent(std::list<int> & calledHaplotypes, 
                                                             std::vector<double> & calledFreqs, 
                                                             std::list<int> candidateHaplotypes, 
                                                             double minLogLik);
        // result
        DindelRealignWindowResult m_result;

        // algorithms

        // initial try: ungapped alignment of read to candidate haplotype, doesn't use base qualities
        void algorithm_hmm(VCFFile & vcfFile);
};

//
// WindowFile
//
class WindowFile
{

    public:
        
        // Constructors
        WindowFile(const std::string & fileName);
        ~WindowFile();

        // Functions
        std::vector<DindelVariant> getNextWindow();
        const std::string & getFilename() const { return m_fileName; }
        int getCurrentLineIndex() const { return m_index; };

    private:

        std::ifstream m_fileHandle;
        bool m_isOpen;
        int m_index;
        std::string m_fileName;
        std::vector<DindelVariant> m_empty;

};

// Add an insertion/deletion into a haplotype
DindelVariant addDeletion(const DindelHaplotype & refHap, const std::string & chrom, 
                          int refPos, int len);
DindelVariant addInsertion(const DindelHaplotype & refHap, const std::string & chrom, 
                           int refPos, const std::string & seq);

//
// CandidateIndelFreq
//
class CandidateIndelFreq
{
    public:

        // Constructors
        CandidateIndelFreq() {};
        CandidateIndelFreq(const std::string & sample, 
                           const std::string & library, 
                           const std::string & fragment) 
        { 
            addObservation(sample, library, fragment); 
        };

        // Functions
        void addObservation(const std::string & sample, const std::string & library, const std::string & fragment);
        int getNumSamples() const { return int(m_samples.size()); };    
        void getStats(int & numLibs, int & numFragments, int & maxFragPerSample) const;

    private:
        class SampleIndels
        {
            public:
                SampleIndels() {};
                void addObservation(const std::string & library, const std::string & fragment);
                typedef HashMap<std::string, int> FragmentToFreq;
                typedef HashMap<std::string, FragmentToFreq > LibToFragment;
                LibToFragment libToFragment;
        };

        HashMap<std::string, SampleIndels> m_samples;
};

typedef std::map<DindelVariant, CandidateIndelFreq> DindelVariantToFreqMap;
typedef std::map<int, DindelVariantToFreqMap > CandidateIndelMap;

#endif  /* DINDELREALIGNWINDOW_H */
