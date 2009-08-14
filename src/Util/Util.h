#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <istream>
#include <cassert>
#include <sstream>

#define CAF_SEP ':'

//
// Enums
//
enum EdgeDir
{
        ED_SENSE = 0,
        ED_ANTISENSE,
        ED_COUNT
};

enum EdgeComp
{
        EC_SAME = 0,
        EC_REVERSE
};

const EdgeDir EDGE_DIRECTIONS[ED_COUNT] = { ED_SENSE, ED_ANTISENSE };

//
// Typedef
//
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::vector<std::string> StringVec;
typedef std::string Sequence;
typedef std::string ContigID;
typedef std::vector<Sequence> SequenceVector;

//
// Alignment
//
struct KAlignment
{
	ContigID contig_id;
	int contig_start_pos;
	int read_start_pos;
	int align_length;
	int read_length;
	bool is_reverse;

	// Return the outer coordinate of the alignment
	int contigOuterCoordinate() const
	{
		if(!is_reverse)
		{
			return contig_start_pos - read_start_pos;
		}
		else
		{
			return contig_start_pos + align_length + read_start_pos;
		}
	}

	//
	// Convert the alignment into the alignment on the reverse complement
	// of the target
	//
	void flipAlignment(int targetLength)
	{
		int tPos = targetLength - contig_start_pos + align_length;
		contig_start_pos = tPos;
		is_reverse = !is_reverse;
	}

	//
	// Get the distance from thto the end of the contig
	// This is in the direction of the alignment
	//
	int getDistanceToEnd(int targetLen) const
	{
		int outerCoordinate = contigOuterCoordinate();
		if(!is_reverse)
			return targetLen - outerCoordinate;
		else
			return outerCoordinate;
	}

	// Comparse by read position
	static int compareReadPos(const KAlignment& a1, const KAlignment& a2)
	{
		return a1.read_start_pos < a2.read_start_pos;
	}

	friend std::istream& operator>> (std::istream& in, KAlignment& a)
	{
		in >> a.contig_id >> a.contig_start_pos;
		in >> a.read_start_pos >> a.align_length;
		in >> a.read_length >> a.is_reverse;
		return in;
	}
};

//
// AlignPair
//
struct AlignPair
{
	friend std::istream& operator>> (std::istream& in, AlignPair& ap)
	{
		std::string readname;
		in >> readname >> ap.aligns[0] >> readname >> ap.aligns[1];
		return in;
	}

	KAlignment aligns[2];
};

//
// AdjInfo
//
struct AdjInfo
{
	ContigID from;
	ContigID to;
	int dir;
	bool comp;

	friend std::istream& operator>>(std::istream& in, AdjInfo& a);
};


//
// Range 
//
struct Range
{
	Range() : start(0), end(0) {}
	Range(int s, int e) : start(s), end(e) {}
	int start;
	int end;

	size_t size() { return end - start; }

	friend std::ostream& operator<<(std::ostream& out, const Range& r)
	{
		out << "[ " << r.start << "," << r.end << " ]";
		return out;
	}

	friend Range intersect(const Range& r1, const Range& r2)
	{
		Range result;
		result.start = std::max(r1.start, r2.start);
		result.end = std::min(r1.end, r2.end);

		// Check for non-overlap
		if(result.end <= result.start)
		{
			result.start = 0;
			result.end = 0;
		}
		return result;
	}
};

//
// Functions
//

//
// Key-value operations
//
template <class C>
std::string makeKeyValue(std::string key, C value)
{
	std::stringstream ss;
	ss << key << CAF_SEP << value;
	return ss.str();
}

StringVec split(std::string in, char delimiter);
void splitKeyValue(std::string in, std::string& key, std::string& value);

//
// Sequence operations
//
Sequence reverseComplement(Sequence seq);
char complement(char base);

//
// Edge Operations
//
inline EdgeDir operator!(const EdgeDir& dir)
{
        return (dir == ED_SENSE) ? ED_ANTISENSE : ED_SENSE;
}

inline EdgeComp operator!(const EdgeComp& comp)
{
        return (comp == EC_SAME) ? EC_REVERSE : EC_SAME;
}

// Correct an edges direction, given the relationship between the nodes
inline EdgeDir correctDir(EdgeDir dir, EdgeComp comp)
{
	return (comp == EC_SAME) ? dir : !dir;
}

#endif

