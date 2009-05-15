#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <istream>
#include <cassert>

//
// Contig description
//
typedef uint32_t ContigID;
struct Contig
{
	ContigID id;
	int length;
	double coverage;
	std::string seq;

	friend std::istream& operator>> (std::istream& in, Contig& c)
	{
		in.ignore(1); // skip ">"
		in >> c.id >> c.length >> c.coverage; // read header
		in.ignore(1); // skip newline
		in >> c.seq; // read seq
		in.ignore(1); // place the ifstream at the next line
		return in;
	}
};

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


typedef std::string Sequence;
Sequence reverseComplement(Sequence seq);
char complement(char base);

#endif

