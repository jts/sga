//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#include <iostream>
#include <math.h>
#include "Util.h"

//
// KAlign
//

// Return the outer coordinate of the alignment
int KAlignment::contigOuterCoordinate() const
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
void KAlignment::flipAlignment(int targetLength)
{
	int tPos = targetLength - contig_start_pos + align_length;
	contig_start_pos = tPos;
	is_reverse = !is_reverse;
}

//
// Get the distance from thto the end of the contig
// This is in the direction of the alignment
//
int KAlignment::getDistanceToEnd(int targetLen) const
{
	int outerCoordinate = contigOuterCoordinate();
	if(!is_reverse)
		return targetLen - outerCoordinate;
	else
		return outerCoordinate;
}

// Comparse by read position
int KAlignment::compareReadPos(const KAlignment& a1, const KAlignment& a2)
{
	return a1.read_start_pos < a2.read_start_pos;
}

// Output
std::istream& operator>> (std::istream& in, KAlignment& a)
{
	in >> a.contig_id >> a.contig_start_pos;
	in >> a.read_start_pos >> a.align_length;
	in >> a.read_length >> a.is_reverse;
	return in;
}

//
// AdjInfo
//

// Input
std::istream& operator>>(std::istream& in, AdjInfo& a)
{
	std::string line;
	getline(in, line);

	// return if we've hit the end
	if(line == "")
		return in;
	
	StringVec fields = split(line, ',');
	assert(fields.size() == 4);

	std::stringstream parser0(fields[0]);
	std::stringstream parser1(fields[1]);
	std::stringstream parser2(fields[2]);
	std::stringstream parser3(fields[3]);

	parser0 >> a.from;
	parser1 >> a.to;
	parser2 >> a.dir;
	parser3 >> a.comp;
	return in;
}

//
// Interval
//
std::ostream& operator<<(std::ostream& out, const Interval& r)
{
	out << r.start << " " << r.end;
	return out;
}

std::istream& operator>>(std::istream& in, Interval& r)
{
	in >> r.start >> r.end;
	return in;
}


Interval intersect(const Interval& r1, const Interval& r2)
{
	Interval result;
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

//
// SeqCoord
//

// Checks if the left edge is extreme
bool SeqCoord::isLeftExtreme() const
{
	return (interval.start == 0 || interval.end == 0);
}

// Checks if the right edge is extreme
bool SeqCoord::isRightExtreme() const
{
	return (interval.end + 1 == seqlen || interval.start + 1 == seqlen);
}

// Checks if at least one endpoint is the end of the string
bool SeqCoord::isExtreme() const
{
	return (isLeftExtreme() || isRightExtreme());
}

// Checks if both end points are the end of the string
bool SeqCoord::isContained() const
{
	return (isLeftExtreme() && isRightExtreme());
}

// Check if the interval is in reverse orientation
bool SeqCoord::isReverse() const
{
	return interval.end < interval.start;
}

// Flip the orientation of the seq coord
void SeqCoord::flip()
{
	//int tmp = interval.end;
	interval.start = seqlen - interval.start - 1;
	interval.end = seqlen - interval.end - 1;
}

std::string SeqCoord::getSubstring(std::string str) const
{
	int left;
	int size; 
	if(interval.start < interval.end)
	{
		left = interval.start;
		size = interval.end - interval.start + 1;
	}
	else
	{
		left = interval.end;
		size = interval.start - interval.end + 1;
	}
	return str.substr(left, size);
}

// Output
std::ostream& operator<<(std::ostream& out, const SeqCoord& sc)
{
	out << sc.id << " " << sc.interval << " " << sc.seqlen;
	return out;
}

// Input
std::istream& operator>>(std::istream& in, SeqCoord& sc)
{
	in >> sc.id >> sc.interval >> sc.seqlen;
	return in;
}

//
// Overlap
//
Overlap::Overlap(std::string i1, int s1, int e1, int l1,
                 std::string i2, int s2, int e2, int l2)
{
	read[0] = SeqCoord(i1, s1, e1, l1);
	read[1] = SeqCoord(i2, s2, e2, l2);
}

// Output
std::ostream& operator<<(std::ostream& out, const Overlap& o)
{
	out << o.read[0] << "\t";
	out << o.read[1];
	return out;
}

// Input
std::istream& operator>>(std::istream& in, Overlap& o)
{
	in >> o.read[0];
	in >> o.read[1];
	return in;
}



//
// Sequence operations
//

// Reverse complement a sequence
Sequence reverseComplement(const Sequence& seq)
{
	std::string out(seq.length(), 'A');
	size_t last_pos = seq.length() - 1;
	for(int i = last_pos; i >= 0; --i)
	{
		out[last_pos - i] = complement(seq[i]);
	}
	return out;
}

// Reverse a sequence
Sequence reverse(const Sequence& seq)
{
	std::string out(seq.length(), 'A');
	size_t last_pos = seq.length() - 1;
	for(int i = last_pos; i >= 0; --i)
	{
		out[last_pos - i] = seq[i];
	}
	return out;
}

// Complement a sequence
Sequence complement(const Sequence& seq)
{
	std::string out(seq.length(), 'A');
	size_t l = seq.length();
	for(size_t i = 0; i < l; ++i)
	{
		out[i] = complement(seq[i]);
	}
	return out;
}

// Complement a base
char complement(char base)
{
	switch(base)
	{
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default:
			assert(false && "Unknown base!");
			return 'N';
	}
}

// Strip the leading directories and
// the last trailling suffix from a filename
std::string stripFilename(std::string filename)
{
	std::string temp(basename(filename.c_str())); // strip leading directory
	size_t suffixPos = temp.find_last_of('.');
	if(suffixPos == std::string::npos)
	{
		return temp; // no suffix
	}
	else
	{
		return temp.substr(0, suffixPos);
	}
}

// Ensure a filehandle is open
void checkFileHandle(std::ifstream& fh, std::string fn)
{
	if(!fh.is_open())
	{
		std::cerr << "Error: could not open " << fn << " for read\n";
		exit(1);
	}	
}


// Split a string into parts based on the delimiter
StringVec split(std::string in, char delimiter)
{
	StringVec out;
	size_t lastPos = 0;
	size_t pos = in.find_first_of(delimiter);

	while(pos != std::string::npos)
	{
		out.push_back(in.substr(lastPos, pos - lastPos));
		lastPos = pos + 1;
		pos = in.find_first_of(delimiter, lastPos);
	}
	out.push_back(in.substr(lastPos));
	return out;
}

// Split a key-value pair
void splitKeyValue(std::string in, std::string& key, std::string& value)
{
	StringVec parts = split(in, CAF_SEP);
	if(parts.size() != 2)
	{
		std::cerr << "Invalid key-value pair " << in << std::endl;
		assert(false);
	}

	key = parts[0];
	value = parts[1];

	assert(key.size() > 0 && value.size() > 0 && "Invalid key-value pair");
}

