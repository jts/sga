#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <istream>
#include <cassert>
#include <sstream>

#define CAF_SEP ':'

//
// Typedef
//
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::string Sequence;
typedef std::string ContigID;

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

	friend std::istream& operator>> (std::istream& in, KAlignment& a)
	{
		in >> a.contig_id >> a.contig_start_pos;
		in >> a.read_start_pos >> a.align_length;
		in >> a.read_length >> a.is_reverse;
		return in;
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

//
// Sequence operations
//
Sequence reverseComplement(Sequence seq);
char complement(char base);


//
// Probability
//
double poisson(int k, double m);
int factorial(int k);

double log_poisson(int k, double m);
double log_factorial(int k);

#endif

