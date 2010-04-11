//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ASQG - Definitions and functions
// for handling ASGQ files
//
#include "ASQG.h"
#include <iterator>

namespace ASQG
{

static char HeaderTag[] = "HT";
static char VertexTag[] = "VT";
static char EdgeTag[] = "ED";

//
RecordType getRecordType(const std::string& record)
{
	// the record type is the first two characters of the record
	if(record.size() < 2)
	{
		std::cerr << "Error: record does not have a valid record-type tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	char recordTag[2];
	record.copy(recordTag, 2);

	if(strncmp(recordTag, HeaderTag, 2) == 0)
	{
		return RT_HEADER;
	}

	if(strncmp(recordTag, VertexTag, 2) == 0)
	{
		return RT_VERTEX;
	}

	if(strncmp(recordTag, EdgeTag, 2) == 0)
	{
		return RT_EDGE;
	}

	// Unknown tag
	std::cerr << "Error: Unknown ASQG file record tag: " << recordTag << std::endl;
	exit(EXIT_FAILURE);
}

//
HeaderRecord parseHeaderRecord(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid header record
	if(tokens[0].compare(HeaderTag) != 0)
	{
		std::cerr << "Error: Record does not have a header tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	HeaderRecord ret;

	// Parse all the tagValues
	std::transform(tokens.begin() + 1, tokens.end(), ret.tags.begin(), SQG::parseTagString);
	return ret;
}

//
VertexRecord parseVertexRecord(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid vertex record
	if(tokens[0].compare(VertexTag) != 0)
	{
		std::cerr << "Error: Record does not have a vertex tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	VertexRecord ret;

	// Parse the ID and sequence
	ret.id = tokens[1];
	ret.seq = tokens[2];

	// Parse all the tagValues
	std::transform(tokens.begin() + 3, tokens.end(), ret.tags.begin(), SQG::parseTagString);
	return ret;
}

//
EdgeRecord parseEdgeRecord(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid vertex record
	if(tokens[0].compare(EdgeTag) != 0)
	{
		std::cerr << "Error: Record does not have an edge tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	EdgeRecord ret;

	// Parse the overlap
	std::stringstream ssparser(tokens[1]);
	ssparser >> ret.overlap;

	// Parse all the tagValues
	std::transform(tokens.begin() + 1, tokens.end(), ret.tags.begin(), SQG::parseTagString);
	return ret;
}

void writeHeaderRecord(std::ostream& out, const HeaderRecord& record)
{
	StringVector fields;
	fields.push_back(HeaderTag);

	std::back_insert_iterator<StringVector> ii(fields);
	std::transform(record.tags.begin(), record.tags.end(), ii, SQG::makeTagString);

	writeFields(out, fields);
}


void writeVertexRecord(std::ostream& out, const VertexRecord& record)
{
	StringVector fields;

	fields.push_back(VertexTag);
	fields.push_back(record.id);
	fields.push_back(record.seq);

	std::back_insert_iterator<StringVector> ii(fields);
	std::transform(record.tags.begin(), record.tags.end(), ii, SQG::makeTagString);
	
	writeFields(out, fields);
}

void writeEdgeRecord(std::ostream& out, const EdgeRecord& record)
{
	StringVector fields;

	fields.push_back(EdgeTag);
	std::stringstream ss;
	ss << record.overlap;
	fields.push_back(ss.str());

	std::back_insert_iterator<StringVector> ii(fields);
	std::transform(record.tags.begin(), record.tags.end(), ii, SQG::makeTagString);
	
	writeFields(out, fields);
}

void writeFields(std::ostream& out, const StringVector& fields)
{
	char delim[2];
	delim[0] = SQG::FIELD_SEP;
	delim[1] = '\0';
	std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, delim));
	out << "\n";
}

}

