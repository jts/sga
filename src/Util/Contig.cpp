#include "Contig.h"

//
// Statics
//
static const char* IDFIELD = "ID";
static const char* LENFIELD = "LN";
static const char* COVFIELD = "CV";
static const char* SEQFIELD = "SQ";
static const char* UNQFIELD = "UN";

Contig::Contig() :	m_length(0), 
					m_coverage(0.0f), 
					m_uniqueFlag(UF_UNKNOWN) {}

std::istream& readFasta(std::istream& in, Contig& c)
{
	in.ignore(1); // skip ">"
	in >> c.m_id >> c.m_length >> c.m_coverage; // read header
	in.ignore(1); // skip newline
	in >> c.m_seq; // read seq
	in.ignore(1); // place the ifstream at the next line
	return in;
}

std::ostream& writeCAF(std::ostream& out, Contig& c)
{
	std::vector<std::string> fields;
		
	// ID
	fields.push_back(makeKeyValue(IDFIELD, c.m_id));
	fields.push_back(makeKeyValue(LENFIELD, c.m_length));
	fields.push_back(makeKeyValue(COVFIELD, c.m_coverage));
	fields.push_back(makeKeyValue(SEQFIELD, c.m_seq));
	fields.push_back(makeKeyValue(UNQFIELD, c.m_uniqueFlag));
	std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, "\t"));
	return out;
}


