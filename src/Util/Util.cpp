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

// Complement the SeqCoord
SeqCoord SeqCoord::complement() const
{
	assert(isExtreme());

	SeqCoord out;
	out.seqlen = seqlen;

	if(isLeftExtreme())
	{
		out.interval.start = std::max(interval.start, interval.end) + 1;
		out.interval.end = out.seqlen - 1;
	}
	else
	{
		out.interval.start = 0;
		out.interval.end = std::min(interval.start, interval.end) - 1;
	}
	assert(out.interval.start <= out.interval.end);
	return out;
}

std::string SeqCoord::getSubstring(const std::string& str) const
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

std::string SeqCoord::getComplementString(const std::string& str) const
{
	SeqCoord comp = complement();
	return comp.getSubstring(str);
}


// Output
std::ostream& operator<<(std::ostream& out, const SeqCoord& sc)
{
	out << sc.interval << " " << sc.seqlen;
	return out;
}

// Input
std::istream& operator>>(std::istream& in, SeqCoord& sc)
{
	in >> sc.interval >> sc.seqlen;
	return in;
}

//
// Sequence operations
//

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
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
std::string reverse(const std::string& seq)
{
	return std::string(seq.rbegin(), seq.rend());
}

// Complement a sequence
std::string complement(const std::string& seq)
{
	std::string out(seq.length(), 'A');
	size_t l = seq.length();
	for(size_t i = 0; i < l; ++i)
	{
		out[i] = complement(seq[i]);
	}
	return out;
}


// count the differences between s1 and s2 over the first n bases
int countDifferences(const std::string& s1, const std::string& s2, size_t n)
{
	int numDiff = 0;
	for(size_t i = 0; i < n; ++i)
	{
		if(s1[i] != s2[i])
			numDiff++;
	}
	return numDiff;
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

// Get the ID of the pair of a given read
std::string getPairID(const std::string& id)
{
	assert(!id.empty());
	std::string pid(id);

	size_t li = id.length() - 1;
	char last = id[li];
	if(last == 'A')
		pid[li] = 'B';
	else if(last == 'B')
		pid[li] = 'A';
	else
	{
		std::cerr << "Unrecognized paired end read id format: " << id << "\n";
		assert(false);
	}
	return pid;
}

// Debug function to parse the distance between two reads
// based on their names
// This assumes is that the read names are just the positions the reads
// were sampled from
size_t debug_getReadDistFromNames(const std::string& name1, const std::string& name2)
{
	std::string id;
	std::string pos;
	StringVec parts = split(name1, ':');
	if(parts.size() != 2)
		return 0;

	int p1 = atoi(parts[1].c_str());
	
	parts = split(name2, ':');
	if(parts.size() != 2)
		return 0;

	int p2 = atoi(parts[1].c_str());
	int dist = int(abs(p1 - p2));
	return dist;
}

