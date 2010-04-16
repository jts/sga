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
		public:
			HeaderRecord();
			
			void setVersionTag(int version);
			void setOverlapTag(int overlapLen);
			void setInputFileTag(const std::string& name);
			void setErrorRateTag(float errorRate);

			void write(std::ostream& out);
			void parse(const std::string& record);

		private:
			
			SQG::IntTag m_versionTag;
			SQG::FloatTag m_errorRateTag;
			SQG::StringTag m_infileTag;
			SQG::IntTag m_overlapTag;
	};

	// A vertex record is an id, sequence and an array of
	// tag:value 
	struct VertexRecord
	{
		public:
			VertexRecord() {}
			VertexRecord(const std::string& i, const std::string& s) : m_id(i), m_seq(s) {}

			void setSubstringTag(bool b);

			void write(std::ostream& out);
			void parse(const std::string& record);

		private:

			std::string m_id;
			std::string m_seq;
			SQG::IntTag m_substringTag;
	};

	// An edge record is just an overlap object and tag:values
	struct EdgeRecord
	{
		public:
			EdgeRecord() {}
			EdgeRecord(const Overlap& o) : m_overlap(o) {}

			void write(std::ostream& out);
			void parse(const std::string& record);

		private:
		
			Overlap m_overlap;
	};

	// Parsing functions
	RecordType getRecordType(const std::string& record);
	HeaderRecord parseHeaderRecord(const std::string& record);
	VertexRecord parseVertexRecord(const std::string& record);
	EdgeRecord parseEdgeRecord(const std::string& record);

	// Writing functions
	void writeFields(std::ostream& out, const StringVector& fields);
};

#endif
