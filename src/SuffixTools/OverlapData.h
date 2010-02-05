//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapData - Data structures for holding overlaps
// between sequence reads
// 
#ifndef OVERLAPDATA_H
#define OVERLAPDATA_H

#include "BWTInterval.h"

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
		bool isTargetRev() const { return data.test(TARGETREV_BIT);	}
		bool isQueryComp() const { return data.test(QUERYCOMP_BIT); }	


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
	OverlapBlock(BWTIntervalPair r, int ol, const BWT* pRB, const AlignFlags& af) : ranges(r), 
	                                                                         overlapLen(ol), 
																			 pRevBWT(pRB), 
																			 flags(af) {}

	// Return the spectrum of extensions given by the interval in ranges
	// The counts are given in the canonical frame, which means that
	// if the query string was reversed, we flip the counts
	AlphaCount getCanonicalExtCount() const;

	static bool sortSizeDescending(const OverlapBlock& ob1, const OverlapBlock& ob2)
	{
		return ob1.overlapLen > ob2.overlapLen;
	}

	BWTIntervalPair ranges;
	int overlapLen;
	const BWT* pRevBWT;
	AlignFlags flags;
};

// The structure written to a file describing the result of the FM-index search for a read
struct OverlapBlockRecord
{
	OverlapBlockRecord() {}
	OverlapBlockRecord(const OverlapBlock& ob) : range(ob.ranges.interval[0]), 
	                                             overlapLen(ob.overlapLen), 
											     flags(ob.flags), 
											     numDiff(0) {}
	// Functions
	friend std::ostream& operator<<(std::ostream& out, const OverlapBlockRecord& obl)
	{
		out << obl.range << " " << obl.overlapLen << " " << obl.flags << " " << obl.numDiff;
		return out;
	}

	friend std::istream& operator>>(std::istream& in, OverlapBlockRecord& obl)
	{
		in >> obl.range >> obl.overlapLen >> obl.flags >> obl.numDiff;
		return in;
	}

	void write(std::ofstream& out)
	{
		range.write(out);
		out.write((char*)&overlapLen, sizeof(overlapLen));
		flags.write(out);
		out.write((char*)&numDiff, sizeof(numDiff));
	}

	void read(std::ifstream& in)
	{
		range.read(in);
		in.read((char*)&overlapLen, sizeof(overlapLen));
		flags.read(in);
		in.read((char*)&numDiff, sizeof(numDiff));
	}

	// Data
	BWTInterval range;
	int overlapLen;
	AlignFlags flags;
	int numDiff;
};

#endif
