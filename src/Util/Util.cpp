//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#include <iostream>
#include <math.h>
#include "Util.h"

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

// 
std::string prefix(const std::string& seq, const unsigned int len)
{
	assert(seq.length() >= len);
	return seq.substr(0, len);
}

// 
std::string suffix(const std::string& seq, const unsigned int len)
{
	assert(seq.length() >= len);
	return seq.substr(seq.length() - len);
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

// count the differences between s1 and s2 over the first n bases
std::string getDiffString(const std::string& s1, const std::string& s2)
{
	std::string out;
	size_t stop = std::min(s1.size(), s2.size());
	for(size_t i = 0; i < stop; ++i)
		out.push_back(s1[i] == s2[i] ? '.' : s1[i]);
	return out;
}

// Strip the leading directories and
// the last trailling suffix from a filename
std::string stripFilename(const std::string& filename)
{
	size_t lastDirPos = filename.find_last_of('/');
	size_t suffixPos = filename.find_last_of('.');
	
	if(lastDirPos == std::string::npos)
		lastDirPos = 0;
	else
		lastDirPos += 1;

	if(suffixPos == std::string::npos)
	{
		return filename.substr(lastDirPos); // no suffix
	}
	else
	{
		return filename.substr(lastDirPos, suffixPos - lastDirPos);
	}
}

// Return the file extension
std::string getFileExtension(const std::string& filename)
{
	size_t suffixPos = filename.find_last_of('.');
	
	if(suffixPos == std::string::npos)
		return "";
	else
		return filename.substr(suffixPos + 1);
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

// Ensure a filehandle is open
void checkFileHandle(std::ofstream& fh, std::string fn)
{
	if(!fh.is_open())
	{
		std::cerr << "Error: could not open " << fn << " for write\n";
		exit(1);
	}	
}


// Split a string into parts based on the delimiter
StringVector split(std::string in, char delimiter)
{
	StringVector out;
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
	StringVector parts = split(in, CAF_SEP);
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
	StringVector parts = split(name1, ':');
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

