#include <iostream>
#include <math.h>
#include "Util.h"

//
// Globals
//

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
	assert(parts.size() == 2 && "Invalid key-value pair");

	key = parts[0];
	value = parts[1];

	assert(key.size() > 0 && value.size() > 0 && "Invalid key-value pair");
}

