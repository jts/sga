//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// EdgeDesc - A unique description of an edge 
//
#include "EdgeDesc.h"

// Operators
bool EdgeDesc::operator<(const EdgeDesc& obj) const
{
	if(id < obj.id)
		return true;
	else if(id > obj.id)
		return false;
	else if(dir < obj.dir)
		return true;
	else if(dir > obj.dir)
		return false;
	else if(comp < obj.comp)
		return true;
	else if(comp > obj.comp)
		return false;
	return false;
}

bool EdgeDesc::operator==(const EdgeDesc& obj) const
{
	return id == obj.id && dir == obj.dir && comp == obj.comp;
}

std::ostream& operator<<(std::ostream& out, const EdgeDesc& ed)
{
	out << ed.id << "," << ed.dir << "," << ed.comp;
	return out;
}

