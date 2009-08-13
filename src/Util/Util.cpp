#include <iostream>
#include <math.h>
#include "Util.h"

//
// Globals
//

//
// AdjInfo
// 
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

