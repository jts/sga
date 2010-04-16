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

const int HEADER_VERSION = 1;

// Record ID tags
static char HEADER_TAG[] = "HT";
static char VERTEX_TAG[] = "VT";
static char EDGE_TAG[] = "ED";

// Header tags
static char VERSION_TAG[] = "VN";
static char OVERLAP_TAG[] = "OL";
static char INPUT_FILE_TAG[] = "IN";
static char ERROR_RATE_TAG[] = "ER";

// Vertex tags
static char SUBSTRING_TAG[] = "SS";

//
// Header Record
//
HeaderRecord::HeaderRecord()
{
	setVersionTag(HEADER_VERSION);
}

//
void HeaderRecord::setVersionTag(int version)
{
	m_versionTag.set(version);
}

//
void HeaderRecord::setOverlapTag(int overlapLen)
{
	m_overlapTag.set(overlapLen);
}

//
void HeaderRecord::setInputFileTag(const std::string& name)
{
	m_infileTag.set(name);
}

//
void HeaderRecord::setErrorRateTag(float errorRate)
{
	m_errorRateTag.set(errorRate);
}

//
void HeaderRecord::write(std::ostream& out)
{
	StringVector fields;
	fields.push_back(HEADER_TAG);

	// Check for mandatory tags
	if(!m_versionTag.isSet())
	{
		std::cerr << "Error: Header version tag not set, aborting." << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		fields.push_back(m_versionTag.toTagString(VERSION_TAG));
	}

	if(m_errorRateTag.isSet())
		fields.push_back(m_errorRateTag.toTagString(ERROR_RATE_TAG));

	if(m_overlapTag.isSet())
		fields.push_back(m_overlapTag.toTagString(OVERLAP_TAG));
	
	if(m_infileTag.isSet())
		fields.push_back(m_infileTag.toTagString(INPUT_FILE_TAG));

	writeFields(out, fields);
}

//
void HeaderRecord::parse(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid header record
	if(tokens[0].compare(HEADER_TAG) != 0)
	{
		std::cerr << "Error: Record does not have a header tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	assert(false);
}

//
// Vertex Record
//
void VertexRecord::setSubstringTag(bool b)
{
	m_substringTag.set(b);
}

//
void VertexRecord::write(std::ostream& out)
{
	StringVector fields;
	fields.push_back(VERTEX_TAG);
	fields.push_back(m_id);
	fields.push_back(m_seq);

	if(m_substringTag.isSet())
		fields.push_back(m_substringTag.toTagString(SUBSTRING_TAG));

	writeFields(out, fields);
}


//
void VertexRecord::parse(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid vertex record
	if(tokens[0].compare(VERTEX_TAG) != 0)
	{
		std::cerr << "Error: Record does not have a vertex tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	m_id = tokens[1];
	m_seq = tokens[2];
	assert(false);
}



//
// EdgeRecord
//
void EdgeRecord::write(std::ostream& out)
{
	StringVector fields;
	fields.push_back(EDGE_TAG);
	std::stringstream ss;
	ss << m_overlap;
	fields.push_back(ss.str());
	writeFields(out, fields);
}

//
void EdgeRecord::parse(const std::string& record)
{
	if(record.size() < 2)
	{
		std::cerr << "Error: Record is not valid\n";
		exit(EXIT_FAILURE);
	}

	// Tokenize record
	StringVector tokens = SQG::tokenizeRecord(record);

	// Ensure the first token indicates this is a valid edge record
	if(tokens[0].compare(EDGE_TAG) != 0)
	{
		std::cerr << "Error: Record does not have an edge tag" << std::endl;
		std::cerr << "Record: " << record << std::endl;
		exit(EXIT_FAILURE);
	}

	// Parse the overlap
	std::stringstream ssparser(tokens[1]);
	ssparser >> m_overlap;
}

//
//
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

	if(strncmp(recordTag, HEADER_TAG, 2) == 0)
	{
		return RT_HEADER;
	}

	if(strncmp(recordTag, VERTEX_TAG, 2) == 0)
	{
		return RT_VERTEX;
	}

	if(strncmp(recordTag, EDGE_TAG, 2) == 0)
	{
		return RT_EDGE;
	}

	// Unknown tag
	std::cerr << "Error: Unknown ASQG file record tag: " << recordTag << std::endl;
	exit(EXIT_FAILURE);
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

