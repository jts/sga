#include "STCommon.h"

//
// Make all the cyclic rotations of a string
//
void makeRotations(SAStringVector& table, SAString s)
{
	int l = s.str.length();
	for(int i = 0; i < l; ++i)
	{
		SAString r = SAString(s.id.label, i, s.str.substr(i, l - i) + s.str.substr(0, i));
		table.push_back(r);
	}
}

