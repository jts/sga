//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SQG - Common definitions and functions for parsing 
// assembly and sequence graph files
//
#ifndef SQG_H
#define SQG_H

#include <vector>
#include <string>
#include "Util.h"

namespace SQG
{
	//
	const char FIELD_SEP = '\t';
	const char TAG_SEP = ':';

	// The types of tags that are possible. 
	// Same as SAM format
	enum TagType
	{
		TT_CHAR = 0,
		TT_INT,
		TT_FLOAT,
		TT_STRING,
		TT_HEX,
		TT_NUM_TYPES
	};

	// Ordering matches above
	const char TagCodes[] = {'A', 'i', 'f', 'Z', 'H'};

	//
	struct TagValue
	{
		static const int NAME_LEN = 3;
		char name[NAME_LEN];
		TagType type;
		std::string value;

		static TagValue makeIntTag(const char* tag, int value);
		static TagValue makeCharTag(const char* tag, char value);
		static TagValue makeStringTag(const char* tag, std::string value);
		static TagValue makeFloatTag(const char* tag, float value);
	};

	typedef std::vector<TagValue> TagValueVector;

	StringVector tokenizeRecord(const std::string& record);
	TagValue parseTagString(const std::string& tagString);
	std::string makeTagString(const TagValue& tv);
};

#endif
