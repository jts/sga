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
// Range
//
std::ostream& operator<<(std::ostream& out, const Range& r)
{
	out << "[ " << r.start << "," << r.end << " ]";
	return out;
}

Range intersect(const Range& r1, const Range& r2)
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

//
// Sequence operations
//

//
// Reverse complement a sequence
//
Sequence reverseComplement(Sequence seq)
{
	std::string out(seq.length(), 'A');
	size_t last_pos = seq.length() - 1;
	for(int i = last_pos; i >= 0; --i)
	{
		out[last_pos - i] = complement(seq[i]);
	}
	return out;
}

//
// Complement a base
//
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
	}
}

//
// Split a string into parts based on the delimiter
//
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

//
// Split a key-value pair
//
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

