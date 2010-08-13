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
#ifndef SCAFFOLDRECORD_H
#define SCAFFOLDRECORD_H

#include "ScaffoldLink.h"
#include "SGUtil.h"

class ScaffoldRecord
{
    public:
        ScaffoldRecord();

        void setRoot(const std::string& root);
        void addLink(const ScaffoldLink& link);

        std::string generateString(const StringGraph* pGraph, bool bNoOverlap, int minOverlap, int maxOverlap, double maxErrorRate) const;

        void parse(const std::string& text);
        void writeScaf(std::ostream* pWriter);

    private:
        
        typedef std::vector<ScaffoldLink> LinkVector;
        
        std::string m_rootID;
        LinkVector m_links;
};

#endif
