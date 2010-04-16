//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ASQG - Definitions and functions
// for handling ASGQ files
//
#ifndef ASQG_H
#define ASQG_H

#include "SQG.h"
#include "Match.h"

namespace ASQG
{
	enum RecordType
	{
		RT_HEADER = 0,
		RT_VERTEX,
		RT_EDGE
	};

	// A header record is just a tag:value pairs
	struct HeaderRecord
	{
		HeaderRecord();
		void addTag(const SQG::TagValue& tag) { tags.push_back(tag); }
		
		void addVersionTag(int version);
		void addOverlapTag(int overlapLen);
		void addInputFileTag(const std::string& name);
		void addErrorRateTag(float errorRate);

		//
		SQG::TagValueVector tags;
	};

	// A vertex record is an id, sequence and an array of
	// tag:value 
	struct VertexRecord
	{
		VertexRecord() {}
		VertexRecord(const std::string& i, const std::string& s) : id(i), seq(s) {}

		//
		void addTag(const SQG::TagValue& tag) { tags.push_back(tag); }
		void addSubstringTag(bool b);

		//
		std::string id;
		std::string seq;
		SQG::TagValueVector tags;
	};

	// An edge record is just an overlap object and tag:values
	struct EdgeRecord
	{
		EdgeRecord() {}
		EdgeRecord(const Overlap& o) : overlap(o) {}
		void addTag(const SQG::TagValue& tag) { tags.push_back(tag); }

		//
		Overlap overlap;
		SQG::TagValueVector tags;
	};


	// Parsing functions
	RecordType getRecordType(const std::string& record);
	HeaderRecord parseHeaderRecord(const std::string& record);
	VertexRecord parseVertexRecord(const std::string& record);
	EdgeRecord parseEdgeRecord(const std::string& record);

	// Writing functions
	void writeHeaderRecord(std::ostream& out, const HeaderRecord& record);
	void writeVertexRecord(std::ostream& out, const VertexRecord& record);
	void writeEdgeRecord(std::ostream& out, const EdgeRecord& record);
	void writeFields(std::ostream& out, const StringVector& fields);
};

#endif
