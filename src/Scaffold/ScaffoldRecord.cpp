//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldRecord - A scaffold consisting of a 
// starting component and a vector of ordered
// links
//
#include "ScaffoldRecord.h"

//
ScaffoldRecord::ScaffoldRecord()
{

}

//
void ScaffoldRecord::setRoot(const std::string& root)
{
    m_rootID = root;
}

//
void ScaffoldRecord::addLink(const ScaffoldLink& link)
{
    m_links.push_back(link);
}

//
void ScaffoldRecord::parse(const std::string& text)
{
    StringVector fields = split(text, '\t');
    assert(fields.size() >= 1);

    m_rootID = fields[0];
    for(size_t i = 1; i < fields.size(); ++i)
    {
        ScaffoldLink link;
        link.parse(fields[i]);
        m_links.push_back(link);
    }
}

//
void ScaffoldRecord::writeScaf(std::ostream* pWriter)
{
    *pWriter << m_rootID;
    for(size_t i = 0; i < m_links.size(); ++i)
        *pWriter << "\t" << m_links[i];
    *pWriter << "\n";
}
