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

	// getTypeCode must be defined for every type that TagValue can take
	static inline char getTypeCode(int)                 { return 'i'; }
	static inline char getTypeCode(char)                { return 'A'; }
	static inline char getTypeCode(const std::string&)  { return 'Z'; }
	static inline char getTypeCode(float)               { return 'f'; }

	// Two-state TagValue that allows the data value to be not set
	template<typename T>
	class TagValue2
	{
		public:
			TagValue2() : m_isSet(false) {}
			TagValue2(const T& v) : m_isSet(true), m_value(v) {}

			T get()
			{
				assert(m_isSet);
				return m_value;
			}

			void set(T v)
			{
				m_isSet = true;
				m_value = v;
			}

			bool isSet()
			{
				return m_isSet;
			}

			std::string toTagString(const char* tag)
			{
				std::stringstream ss;
				char type_code = getTypeCode(m_value);
				ss << tag << TAG_SEP << type_code << TAG_SEP << m_value;
				return ss.str();
			}
		
		private:
			T m_value;
			bool m_isSet;
	};

	// These are the valid tags that can be used
	typedef TagValue2<char> CharTag;
	typedef TagValue2<int> IntTag;
	typedef TagValue2<float> FloatTag;
	typedef TagValue2<std::string> StringTag;

/*
	class CharTag
	{
		public:
			bool isSet();
			int get();
			void set(int v);
			std::string toString();

		private:
			int value;
			bool m_isSet;
	}	

	class IntTag
	{
		public:
			bool isSet();
			int get();
			void set(int v);
			std::string toString();

		private:
			int value;
			bool m_isSet;
	}

	class FloatTag
	{
		public:
			bool isSet();
			void set(float v);
			float get();
			std::string toString();

		private:
			float value;
			bool m_isSet;
	}

	class StringTag
	{
		public:
			bool isSet();
			void set(std::string v);
			float get();
			std::string toString();

		private:
			std::string value;
			bool m_isSet;
	}
*/

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
