//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SequenceWorkItem - Definition of the data structure used in the generic
// functions to process a sequence
//
#ifndef SEQUENCEWORKITEM_H
#define SEQUENCEWORKITEM_H

struct SequenceWorkItem
{
    SequenceWorkItem(size_t ri, const SeqRecord& sr) : idx(ri), read(sr) {}
    size_t idx;
    SeqRecord read;
};

#endif
