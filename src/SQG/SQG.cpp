//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SQG - Common functions for parsing 
// assembly and sequence graph files
//
#include "SQG.h"

namespace SQG
{

//
StringVector tokenizeRecord(const std::string& record)
{
	return split(record, FIELD_SEP);
}

//
TagValue TagValue::makeIntTag(const char* tag, int value)
{
	TagValue ret;
	strncpy(ret.name, tag, NAME_LEN);

	ret.type = TT_INT;

	std::stringstream ss;
	ss << value;
	ret.value = ss.str();
	return ret;
}

//
TagValue TagValue::makeCharTag(const char* tag, char value)
{
	TagValue ret;
	strncpy(ret.name, tag, NAME_LEN);

	ret.type = TT_CHAR;

	std::stringstream ss;
	ss << value;
	ret.value = ss.str();
	return ret;
}

//
TagValue TagValue::makeStringTag(const char* tag, std::string value)
{
	TagValue ret;
	strncpy(ret.name, tag, NAME_LEN);

	ret.type = TT_STRING;

	std::stringstream ss;
	ss << value;
	ret.value = ss.str();
	return ret;
}

//
TagValue TagValue::makeFloatTag(const char* tag, float value)
{
	TagValue ret;
	strncpy(ret.name, tag, NAME_LEN);

	ret.type = TT_FLOAT;

	std::stringstream ss;
	ss << value;
	ret.value = ss.str();
	return ret;
}


//
TagValue parseTagString(const std::string& tagString)
{
	StringVector tokens = split(tagString, TAG_SEP);
	if(tokens.size() != 3)
	{
		std::cerr << "Error: Invalid Tag - expected 3 fields" << std::endl;
		std::cerr << "Tag: " << tagString << std::endl;
		exit(EXIT_FAILURE);
	}

	// Copy the tag code
	TagValue ret;
	tokens[0].copy(ret.name, TagValue::NAME_LEN);

	// Parse the type of the tag
	bool bValidType = false;
	if(tokens[1].size() == 1)
	{
		for(size_t i = 0; i < TT_NUM_TYPES; ++i)
		{
			if(tokens[1][0] == TagCodes[i])
			{
				ret.type = static_cast<TagType>(i);
				bValidType = true;
				break;
			}
		}
	}

	if(!bValidType)
	{
		std::cerr << "Error: Invalid Tag type" << std::endl;
		std::cerr << "Tag type: " << tokens[1] << std::endl;
		exit(EXIT_FAILURE);
	}
	
	ret.value = tokens[2];
	return ret;
}

// Convert a tag to a string
std::string makeTagString(const TagValue& tv)
{
	std::stringstream ss;
	ss << tv.name << TAG_SEP << TagCodes[tv.type] << TAG_SEP << tv.value;
	return ss.str();
}

}
