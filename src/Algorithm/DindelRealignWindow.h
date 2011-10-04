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
#define	DINDELREALIGNWINDOW_H
#include <fstream>
#include "Util.h"
#include "HashMap.h"
#include "MultiAlignment.h"
#include <iomanip>
#include <list>
#include <set>
#include <map>
enum DindelWindowCandHapAlgorithm {LINEAR} ;

const size_t DINDEL_HASH_SIZE=8;
const int DINDEL_HMM_BANDWIDTH=6;

typedef std::string SampleName;

long long int combinations(const int n, const int k);

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

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

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


std::vector<std::string> SplitString(const std::string & str, char sep);

void parseRegionString(const std::string & region, std::string & chrom, int & start, int & end);

class DindelRead;



class VariantPriors
{
    public:
        VariantPriors();
        double getDefaultProbVariant(const std::string & type) const;
        double m_probSNP, m_probINDEL, m_probMNP;

};



// wrapper around BAM
class DindelRead
{
    public:
        DindelRead(const SeqItem* pSeqItem, const SampleName & sampleName, double mappingQual, int fixedBaseQual);
        const std::string  getSequence() const { return m_pSeqItem->seq.toString(); }
        const std::vector<double> getBaseQuals() const;
        const std::string getID() const { return m_pSeqItem->id; }
        double getMappingQual() const { return m_mappingQual; } //return double(bam->core.qual); }
        double getLogProbNotMapping() const { return m_mappingQual*0.23026; }
	bool isForward() const { return true; } //FIXME
        bool mateIsForward() const { return false; } //FIXME
        //BAM bool BAMCigarHasIndel() const;
        //BAM int getBamStartPos() const { return bam->core.pos; }
        //BAM int getBamStartPosAdjusted() const;
        //BAM int getBAMEnd() const { return bam->core.n_cigar? bam_calend(&bam->core, bam1_cigar(bam)) : bam->core.pos + 1; }
        //BAM int getBAMEndAdjusted() const;
        bool isUnmapped() const { return false; } //FIXME
        bool mateIsUnmapped() const { return false; } //FIXME
        char getBase(int b) const { return m_pSeqItem->seq.get(size_t(b)); }
        int getQual(int b) const { b=1; return m_fixedBaseQual; }
        int length() const { return int(m_pSeqItem->seq.length()); }
        //BAM const bam1_t *getBam() const { return bam; }

        void setMappingQuality(double mappingQual) { m_mappingQual=mappingQual; }
        void setRCRead(bool status) { m_rcRead=status; if (status) setupHash(); } // need to recompute hash
        bool getRCRead() const { return m_rcRead; }
        void getLogProbCorrectError(std::vector<double> & lpCorrect, std::vector<double> & lpError) const;
        const std::vector<unsigned int> & getHashKeys();
        const std::string & getSampleName() const { return m_sampleName; }
        const std::string getLibraryName() const { return std::string("NAN"); }//FIXME
    private:
        //BAM const bam1_t *bam;
        //BAM const DindelBAM *m_pDindelBAM;
        double m_mappingQual;
        int m_fixedBaseQual;
        bool m_rcRead, m_setupHash; // reverse-complement the read. Intended for unmapped reads. Determined using the status of the mate.
        SampleName m_sampleName;
        const SeqItem *m_pSeqItem;
        
        std::vector<unsigned int> m_hashKeys; // stores hashes for the read

        void setupHash();
};


class DindelVariant
{
    public:
  	DindelVariant();      
        DindelVariant(const std::string & chrom, const std::string & ref, const std::string & alt, int pos);

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

        void setup();
        void checkSequence(const std::string & seq); // checks whether sequence has no funny characters

};

// sequence hash

class DindelSequenceHash
{
    public:
	DindelSequenceHash() {};
	DindelSequenceHash(const std::string & sequence);

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
	typedef std::tr1::unordered_map<unsigned int, std::list<int> > Hash;
        typedef std::tr1::unordered_map<unsigned int, int > UniqueHash;


        Hash m_hash;
        // UniqueHash m_uniqueHash;
	std::set<int> emptySet;
        std::list<int> emptyList;

        void makeHash(const std::string & sequence);
};

const int DELETION_NOVEL_SEQUENCE = -1;
const int SNP = -2;
const int MULTINUCLEOTIDE_RUN = -3;
const int INSERTION = -4;
const int LEFTOVERHANG = -5;
const int RIGHTOVERHANG = -6;
const int ADDVARIANT_DEBUG = -9;

class DindelHaplotype
{
    public:
        // a DindelHaplotype must be constructed starting from the reference haplotype.
        // Differences with the reference can be added by calling addVariant
        DindelHaplotype(const std::string & refSeq, int refSeqStart, bool isReference);
        DindelHaplotype(const std::string & refName, const std::string & refSeq, int refSeqStart, const MultiAlignment & ma, int varRow, int refRow);
	DindelHaplotype(const DindelHaplotype & haplotype, int copyOptions);
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
        // initializes haplotype from reference sequence
        // sets homopolymer length vector
        void initHaplotype();
        void determineHomopolymerLengths();
        void updateHaplotype();
        
        std::string m_seq;
        std::vector<int> m_hplen; // gives the homopolymer run length for any base in the haplotype sequence going left and right.
        std::vector<int> m_refPos; // position of haplotype base on reference.

        // list of VCF4 style variants contained in haplotype
        std::vector<DindelVariant> m_variants;
        std::tr1::unordered_map<std::string, std::pair<int, int> > m_variant_to_pos;

        // position of leftmost haplotype base on reference sequence
        int m_hapRefStart;

	bool m_isReference;
        // constants for m_refPos

        DindelSequenceHash m_sequenceHash;
};





// DindelWindow contains a list of variants, the candidate haplotypes generated from those variants, and the relationships
// between the candidate haplotypes and the variants.
// The candidate haplotypes are annotated with respect to the reference sequence.
// The candidate variants are specified as modifications to the reference sequence by any string.replace operation
class DindelWindow
{
    public:

        DindelWindow(const std::vector<DindelVariant> & variants, const std::string & refHap, int refHapStart, int windowPad);
        //~DindelWindow();

        // Create window from a set of haplotypes and a reference sequence.
        // Uses SGA MultiAlignment to annotate the variations in the haplotypes with respect to the reference sequence.
        DindelWindow(const std::vector<std::string> & haplotypeSequences, const std::string & refHap, int refHapStart, const std::string & refName);
	~DindelWindow();
	const std::vector<DindelHaplotype> & getHaplotypes() const { return haplotypes; }
        const std::string getChrom() const { return m_chrom; }
        int getLeft() const { return m_leftRef; }
        int getRight() const { return m_rightRef; }
        void addVariant(const DindelVariant & variant, bool addToAll);
        void filterHaplotypes(const std::vector<bool> & filterHaplotype);

    private:
        void makeWindow(const std::vector<DindelVariant> & variants);
	DindelHaplotype makeAltHapFromReference(const DindelVariant & var );
        void generateCandidateHaplotypes(const std::vector<DindelVariant> & variants);
        void setupRefSeq(const std::string & refHap, int refHapStart);
	bool addHaplotype(const DindelHaplotype & haplotype);
        void generateCandidateHaplotypesFromMultiAlignment();

	std::vector<DindelHaplotype> haplotypes; // candidate haplotypes
      
        // leftmost and rightmost coordinates of window on reference sequence
        int m_leftRef, m_rightRef;
        std::string m_refSeq; // reference sequence for [m_leftRef, m_rightRef]
        std::string m_chrom;

        // links variant to haplotype
        // std::vector<std::vector<int> > m_variantToHaplotype;

        // keeps track of haplotypes.
        std::map<std::string, int> m_hashAltHaps;

        //Fasta *m_pFasta; // pointer to Fasta class for getting the reference sequence

        //parameters

        int m_windowPad;
        DindelWindowCandHapAlgorithm m_candHapAlgorithm;
        MultiAlignment *m_pHaplotype_ma;
};
class DindelRealignWindow;

class VCFFile
{
    public:
        class VCFEntry
        {
	    friend class VCFFile;
            public:
                typedef std::map<std::string, std::string> Maps;
                VCFEntry() { m_isEmpty = true; };
                std::string chrom, ref, alt, id, infoString;
                int pos, qual;
                Maps  info, filters;
                bool isEmpty() { return m_isEmpty; };
                std::string getInfoValue(const std::string & key) const;
                std::vector<std::string> getFilters() const;
                void parseInfoString();
                template<class T> bool fromInfoTag(T & val, const std::string & tag, std::ios_base& (*f)(std::ios_base&));
		void setFilters(const std::string & str);
		void write(std::ostream & out);
            private:
                bool m_isEmpty;

        };


	 VCFFile();
	 VCFFile(const std::string & fileName, const std::string & mode);
	 std::ostream & getOutputStream() { assert(m_mode=="w"); return m_outputFileHandle; }
         VCFEntry getNextEntry();
         void setSamples(const std::vector<std::string> & samples);
         const std::vector<std::string> & getSamples() const { return m_samples; }
	 ~VCFFile();
         void outputHeader(const std::string & refFile, const std::string & paramString);
    private:
	 std::ofstream m_outputFileHandle;
	 std::ifstream m_inputFileHandle;
	 std::string m_fileName;
	 bool m_isOpen;
         std::string m_mode;
         std::vector<std::string> m_samples;
};

struct CoalescentResult
{
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
    double qual, singleSampleQual;
    double mapFreq; // MAP frequency of haplotype
    bool singleSampleCall; // call based on simply thresholding individual sample genotype/haplotype likelihoods.
};

class DindelRealignWindowResult
{
    public:
        class Inference
        {
            public:
                //  phred-scaled posterior prob
                double qual;
		double freq;
                double strandBias;
                // std::set<int> readsForward, readsForwardZeroMismatch;
		// std::set<int> readsReverse, readsReverseZeroMismatch;
                int numReadsForward, numReadsReverse, numReadsForwardZeroMismatch, numReadsReverseZeroMismatch;
		int numUnmapped;
                int numLibraries, numReadNames;
                
		bool isSingleSampleCall;
		double singleSampleQual;

                std::vector<int> histDistance; // histogram for distances to variant
                std::vector<int> histAlignLik; // histogram for read-haplotype alignment likelihoods
                std::vector<int> histMapQ; // histogram for mapping qualities
                std::string infoStr; // special info string that can be set freely
                void outputAsVCF(const DindelVariant & var, const DindelRealignWindowResult & result, VCFFile & vcfFile) const;
                static double computeStrandBias(int numForward, int numReverse);
                void addDistanceToHistogram(int distance);  // these functions will take care of scaling and initialization
                void addAlignLikToHistogram(double logLik); //
                void addMapQToHistogram(double mappingQuality);
        };
        class GenotypeCall
        {
            public:
                GenotypeCall() { qual=0.0; count = 0; called = false;}
                double qual;
                int count;
                bool called; // variant was present in called pair of haplotypes
                double gl[3];
        };

	typedef std::map<DindelVariant, Inference> VarToInference;
        typedef std::map<DindelVariant, GenotypeCall> VarToGenotypeCall;
        typedef std::tr1::unordered_map<std::string, VarToGenotypeCall> SampleToGenotypes;

   	DindelRealignWindowResult() : m_pDindelRealignWindow(NULL){ };

        DindelRealignWindowResult(const DindelRealignWindow & dindelRealignWindow) : m_pDindelRealignWindow(&dindelRealignWindow){ };
	
	std::vector<DindelHaplotype> haplotypes;

        std::vector<double> haplotypeFrequencies;
        VarToInference variantInference; // integrates same variants across different haplotypes
        SampleToGenotypes sampleToGenotypes;
        
	int numHapSpecificReads; // Number of reads that were used to call variants (that map well to one haplotype but not another)
        int numHaplotypes; // number of inferred haplotypes
        bool outputGenotypes;
        void outputVCF(VCFFile & vcfFile);
        
    private:
    
        const DindelRealignWindow *m_pDindelRealignWindow;
};



class CandidateSNP
{
    public:
	CandidateSNP(char refBase)
	{
	    m_refBase = refBase;
	    baseToReads = BaseToReads(4);
	}
	char getBase(size_t idx) const { return "ACGT"[idx]; }
	size_t getBaseIdx(char base) const { switch(base) { case 'A' : return 0; case 'C' : return 1; case 'G' : return 2; case 'T' : return 3; case 'N' : throw std::string("N not allowed"); }; return 0; }
	void addRead(char refBase, char altBase, size_t readIndex)
	{
	    assert(m_refBase==refBase);
	    baseToReads[getBaseIdx(altBase)].insert(readIndex);
	}
	typedef std::vector < std::set<size_t> > BaseToReads;
	BaseToReads baseToReads;
	char m_refBase;
};


class PosToCandidateSNP : public std::tr1::unordered_map<int, CandidateSNP >
{
    public:
	void addSNP(int readIndex, int refPos, char refBase, char altBase);

};

class ReadHaplotypeAlignment
{
    public:

        ReadHaplotypeAlignment() { logLik = 1.0; nmm =-1; postProbLastReadBase=-1.0; isUngapped=false; hapPosLastReadBase=-1;}
        ReadHaplotypeAlignment(double _logLik, int _nmm) { logLik = _logLik; nmm = _nmm; postProbLastReadBase=-1.0; isUngapped=false; hapPosLastReadBase=-1;}
        ReadHaplotypeAlignment(double _logLik, int _nmm, bool _isUngapped) { logLik = _logLik; nmm = _nmm; postProbLastReadBase=-1.0; isUngapped=_isUngapped; hapPosLastReadBase=-1;}
       
        double logLik;
        double postProbLastReadBase; // posterior probability of last read base
        int nmm; // number of mismatches between read and haplotype
	int hapPosLastReadBase; // position of last readbase on haplotype. -1 indicates off haplotype (shouldn't happen)
        bool isUngapped;
};

class DindelRealignParameters
{
    friend class DindelRealignAlgorithm;
    friend class DindelRealignWindow; 
    private:
       std::string getDefaultParameters() const;
       void initFromString(const std::string & paramString);

    public:
        DindelRealignParameters();
        DindelRealignParameters(const std::string & paramString);// { initFromString(paramString); }
    private:
        const std::string & getParamString() const { return m_paramString; }
        void checkAndInit();

        std::string m_paramString;

        std::string excludeSamplesFile;

        int windowReadBuffer;
        int minVariantSep, haplotypeWidth, minCandidateAlleleCount;
        int maxNumReads; // maximum number for fetch
        int maxNumReadsWindow; // maximum number of reads in window
        int maxNumCandidatesPerWindow;
        double theta_snp, theta_indel;

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
};

// Realigns reads in a given window.
// pOverlapper and dindelReads should correspond 1-1. It is assumed that pOverlapper gives back indices into the dindelReads array.
class DindelRealignWindow
{
    public:
        
        DindelRealignWindow(const DindelWindow *pDindelWindow, std::vector<DindelRead> & dindelReads, const DindelRealignParameters & realignParameters);

        void run(const std::string & algorithm, VCFFile & vcfFile);

        //~DindelRealignWindow();

        

        
	const DindelWindow & getDindelWindow() const { return  m_dindelWindow; }
        const std::vector< std::vector< ReadHaplotypeAlignment > > & getHapReadAlignments() const { return hapReadAlignments; }
    private:
        DindelWindow m_dindelWindow;
        std::vector<DindelRead> *m_pDindelReads;
	//const DindelOverlapper* m_pDindelOverlapper;
        //const Fasta *m_pFasta;
        DindelRealignParameters realignParameters;
        std::vector<int> m_readMapsToReference;

        // realignment helper classes

      
        // HAPLOTYPE ALIGNMENT BUSINESS

        PosToCandidateSNP m_posToCandidateSNP;
        void filterHaplotypes();
        void addSNPsToHaplotypes();
	void processHaplotypes(size_t firstHapIdx, size_t lastHapIdx);
        // convert seed positions into a likelihood given the candidate haplotype
	ReadHaplotypeAlignment computeReadHaplotypeAlignment(size_t readIndex, const DindelHaplotype& haplotype, int start, int end, const std::vector<double> & lpError, const std::vector<double> & lpCorrect, bool rcReadSeq);
	void HMMAlignReadAgainstHaplotypes(size_t readIndex, size_t firstHap, size_t lastHap, const std::vector<double> & lpCorrect, const std::vector<double> & lpError);

        // HAPLOTYPE FREQUENCY ESTIMATION BUSINESS

	// haplotype frequencies
	class HL
   	{
    	    public:
    		HL(double _ll, int _idx) { ll = _ll; idx = _idx;}
    	        bool operator<(const HL & hl) const { 
                    if (int(ll*10.0)==int(hl.ll*10.0))
                    {
                        if (idx<hl.idx) return true; else return false;

                    } else return int(ll*10.0)<int(hl.ll*10.0);
		}
     	        double ll;
    	        int idx;
    	};

        void printReadAlignments(int readIdx, std::ostream & out, int offset, bool supportAlt);
    	DindelRealignWindowResult estimateHaplotypeFrequencies(double minLogLikAlignToRef, double minLogLikAlignToAlt, bool capUsingMappingQuality, bool print);
	double getHaplotypePrior(const DindelHaplotype & h1, const DindelHaplotype & h2) const;
        // this is a function that hides how the read information is actually stored in the various ReadTables
        const DindelRead & getRead(size_t readIndex) const { return (*m_pDindelReads)[readIndex]; }

	// haplotype-read likelihoods
	std::vector< std::vector< ReadHaplotypeAlignment > > hapReadAlignments; // first dim is haplotypes, second dim is reads
	// update hapReadAlignments.

        void computeReadHaplotypeAlignmentsUsingHMM(size_t firstHap, size_t lastHap);
        void addDiploidGenotypes(DindelRealignWindowResult & result, bool useEstimatedHaplotypeFrequencies);

        
        typedef std::tr1::unordered_map<int, CoalescentResult> HapIdxToCoalescentResult;

        class ReadSamples
        {
            public:
                ReadSamples(const std::vector<DindelRead> & dindelReads);
                std::vector<std::string> samples;
                std::vector<int> readIndices;
		std::vector<double> logNumStates;
            private:
                const std::vector<DindelRead> *m_pDindelReads;
                


        };

        ReadSamples m_readSamples;

        HapIdxToCoalescentResult computeHaplotypesCoalescent(std::list<int> & calledHaplotypes, std::vector<double> & calledFreqs, std::list<int> candidateHaplotypes, double minLogLik);
        // result
        DindelRealignWindowResult m_result;

        // algorithms

        // initial try: ungapped alignment of read to candidate haplotype, doesn't use base qualities
        void algorithm_hmm(VCFFile & vcfFile);
};


class WindowFile
{

    public:
        WindowFile(const std::string & fileName);
        std::vector<DindelVariant> getNextWindow();
        const std::string & getFilename() const { return m_fileName; }
        int getCurrentLineIndex() const { return m_index; };
        ~WindowFile();

    private:
        std::ifstream m_fileHandle;
	bool m_isOpen;
        int m_index;
        std::string m_fileName;
        std::vector<DindelVariant> m_empty;

};






/*

 Functions for extracting candidate indels from a set of BAM files
 
 */
DindelVariant addDeletion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, int len);
DindelVariant addInsertion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, const std::string & seq);

class CandidateIndelFreq
{
public:
    CandidateIndelFreq() {  };
    CandidateIndelFreq(const std::string & sample, const std::string & library, const std::string & fragment) { addObservation(sample, library, fragment); };
    void addObservation(const std::string & sample, const std::string & library, const std::string & fragment);
    int getNumSamples() const { return int(m_samples.size()); };    
    void getStats(int & numLibs, int & numFragments, int & maxFragPerSample) const;
private:
    class SampleIndels
    {
	public:
	SampleIndels() {};
        void addObservation(const std::string & library, const std::string & fragment);
        typedef std::tr1::unordered_map<std::string, int> FragmentToFreq;
        typedef std::tr1::unordered_map<std::string, FragmentToFreq > LibToFragment;
        LibToFragment libToFragment;
    };
    std::tr1::unordered_map<std::string, SampleIndels> m_samples;

    
};

typedef std::map<DindelVariant, CandidateIndelFreq> DindelVariantToFreqMap;
typedef std::map<int, DindelVariantToFreqMap > CandidateIndelMap;






#endif	/* DINDELREALIGNWINDOW_H */

