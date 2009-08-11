#ifndef STCOMMON_H
#define STCOMMON_H

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
#include <set>

struct Ident
{
	Ident() : label(""), idx(0) {}
	Ident(std::string l, int i) : label(l), idx(i) {}
	Ident(std::string l) : label(l), idx(0) {}

	friend std::ostream& operator<<(std::ostream& out, const Ident& i)
	{
		out << i.label << "," << i.idx;
		return out;
	}

	std::string label;
	int idx;
};

//
//
//
class SAString
{	
	public:

		SAString(std::string m, int i, std::string s) : id(m,i), str(s) {}
		SAString(std::string m, std::string s) : id(m), str(s) {}
		
		friend int operator<(const SAString& o1, const SAString& o2)
		{
			return o1.str < o2.str;
		}

		friend std::ostream& operator<<(std::ostream& out, const SAString& s)
		{
			out << s.id << "\t" << s.str;
			return out;
		}

		Ident id;
		std::string str;
};

//
// Typedefs
//
typedef std::vector<SAString> SAStringVector;
typedef std::vector<int> IntVector;
typedef std::vector<Ident> IdentVector;
typedef std::map<char, int> CharIntMap;
typedef std::set<char> CharSet;
typedef std::map<std::string, int> IDIntMap;

//
// Make all the cyclic rotations of a string
//
void makeRotations(SAStringVector& table, SAString s);


//
// Print out a map using cout
//
template<class K, class V>
void printMap(const std::map<K,V>& m)
{
	for(typename std::map<K,V>::const_iterator iter = m.begin(); iter != m.end(); ++iter)
	{
		std::cout << iter->first << "\t" << iter->second << "\n";
	}
}

//
// Print a vector
//
template<class T>
void printVector(const std::vector<T>& v)
{
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, "\n"));
}

#endif
