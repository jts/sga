//-------------------------------------------------------------------------------
// 
// MultipleAlignment - Class for constructing and managing a multiple alignment
// constructed from a set of pairwise overlaps
//
// Copyright (C) 2011 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "multiple_alignment.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <assert.h>
#include <stdio.h>

//#define MA_DEBUG 1

//
// MultipleAlignmentElement
//

//
MultipleAlignmentElement::MultipleAlignmentElement(const std::string& _name, 
                                                   const std::string& _sequence,
                                                   const std::string& _quality,
                                                   size_t leading,
                                                   size_t trailing) : name(_name), 
                                                                      padded_sequence(_sequence),
                                                                      padded_quality(_quality),
                                                                      leading_columns(leading),
                                                                      trailing_columns(trailing)
{

}

//
size_t MultipleAlignmentElement::getNumColumns() const
{
    return leading_columns + padded_sequence.size() + trailing_columns;
}

//
char MultipleAlignmentElement::getColumnSymbol(size_t column_idx) const
{
    assert(column_idx < getNumColumns());
    if(column_idx < leading_columns || column_idx >= leading_columns + padded_sequence.size()) {
        return '\0';
    }
    else {
        assert(column_idx - leading_columns < padded_sequence.size());
        return padded_sequence[column_idx - leading_columns];
    }
}

//
char MultipleAlignmentElement::getColumnQuality(size_t column_idx) const
{
    assert(column_idx < getNumColumns());
    if(padded_quality.empty() || 
       column_idx < leading_columns || 
       column_idx >= leading_columns + padded_quality.size()) {
        return '\0';
    }
    else {
        assert(column_idx - leading_columns < padded_sequence.size());
        return padded_quality[column_idx - leading_columns];
    }
}


//
int MultipleAlignmentElement::getPaddedPositionOfBase(size_t idx) const
{
    size_t unpadded_count = 0;
    for(size_t i = 0; i < padded_sequence.size(); ++i) {
        if(padded_sequence[i] != '-') {
            if(unpadded_count == idx)
                return i;
            else
                unpadded_count += 1;
        }
    }
    std::cerr << "Base index out of bounds: " << idx << "\n";
    assert(false);
    return -1;
}

//
void MultipleAlignmentElement::insertGapBeforeColumn(size_t column_index)
{
    // Check if the column to insert the gap falls within the leading columns
    // If so, just increase the number of leading columns to account for the inserted base.
    // If the column index is one greater than the offset, then we want to 
    // insert a gap before the first base. This is equivalent to just
    // extending the offset.
    if(column_index <= leading_columns + 1) {
        leading_columns += 1;
    }
    else {
        assert(column_index > leading_columns);
        size_t insert_position = column_index - leading_columns;
        if(insert_position < padded_sequence.size()) {
            padded_sequence.insert(insert_position, 1, '-');
            if(!padded_quality.empty())
                padded_quality.insert(insert_position, 1, '-');
        }
        else
            trailing_columns += 1;
    }
}

//
std::string MultipleAlignmentElement::getUnpaddedSequence() const
{
    std::string out;
    for(size_t i = 0; i < padded_sequence.size(); ++i) {
        if(padded_sequence[i] != '-')
            out.push_back(padded_sequence[i]);
    }
    return out;
}

//
// MultipleAlignment
//

//
void MultipleAlignment::addBaseSequence(const std::string& name, const std::string& sequence, const std::string& quality)
{
    m_sequences.push_back(MultipleAlignmentElement(name, sequence, quality, 0, 0));
}

// See header
void MultipleAlignment::addOverlap(const std::string& incoming_name,
                                   const std::string& incoming_sequence,
                                   const std::string& incoming_quality,
                                   const SequenceOverlap& reference_incoming_overlap)
{
    // This function cannot be called before a base element has been added
    assert(!m_sequences.empty());
    MultipleAlignmentElement* template_element = &m_sequences.front();
    _addSequence(incoming_name, incoming_sequence, incoming_quality, template_element, reference_incoming_overlap);
}

//
void MultipleAlignment::addExtension(const std::string& incoming_name,
                                     const std::string& incoming_sequence,
                                     const std::string& incoming_quality,
                                     const SequenceOverlap& previous_incoming_overlap)
{
    // This function cannot be called before a base element has been added
    assert(!m_sequences.empty());
    MultipleAlignmentElement* template_element = &m_sequences.back();
    _addSequence(incoming_name, incoming_sequence, incoming_quality, template_element, previous_incoming_overlap);
}

// Adds a new string into the multiple alignment using the overlap
// between the incoming sequence and an existing sequence in the
// multiple alignment to calculate the new padded string
void MultipleAlignment::_addSequence(const std::string& name,
                                     const std::string& sequence,
                                     const std::string& quality, 
                                     MultipleAlignmentElement* template_element, 
                                     const SequenceOverlap& overlap)
{
    // Get the padded sequence for the template element
    const std::string& template_padded = template_element->padded_sequence;

    // The output padded sequence for the incoming
    std::string padded_output;
    std::string padded_quality;
    
    // Sanity checks
    assert(quality.empty() || quality.size() == sequence.size());

    // Iterate over the cigar string and the padded sequence of the template element
    // to determine where to insert gap symbols
    size_t cigar_index = 0;
    size_t template_index = template_element->getPaddedPositionOfBase(overlap.match[0].start);
    size_t incoming_index = overlap.match[1].start;
    size_t template_leading = template_element->leading_columns;
    size_t incoming_leading = template_index + template_leading;

    // Expand the cigar for easier parsing
    std::string expanded_cigar = expandCigar(overlap.cigar);
    assert(!expanded_cigar.empty());
    assert(template_index < template_padded.size());
    assert(template_padded[template_index] != '-');

#ifdef MA_DEBUG
    std::cout << "Cigar: " << expanded_cigar << "\n";
    std::cout << "template: " << template_padded.substr(template_index) << "\n";
    std::cout << "incoming: " << sequence.substr(incoming_index) << "\n";
    std::cout << "Pairwise:\n";
    overlap.printAlignment(template_element->getUnpaddedSequence(), sequence);
#endif

    while(cigar_index < expanded_cigar.size()) {

        // Check if we are in an existing template gap. This must be handled
        // seperately
        bool in_template_gap = template_padded[template_index] == '-';
        if(in_template_gap) {
            // If we are in a incoming sequence insertion
            // (cigar D) then we are adding a base into a known
            // gap. Add the current incoming base to the output
            if(expanded_cigar[cigar_index] == 'D') {
                padded_output.push_back(sequence[incoming_index]);
                if(!quality.empty())
                    padded_quality.push_back(quality[incoming_index]);

                incoming_index += 1;
                cigar_index += 1;
                template_index += 1;
            } else { 
                // This is an insertion that is in some other sequence
                // in the multiple alignment. Add a gap to the padded output
                padded_output.push_back('-');
                if(!quality.empty())
                    padded_quality.push_back('-');
                
                // Increment the template index
                template_index += 1;    
            }
        } else {
            // Not a template gap
            switch(expanded_cigar[cigar_index]) {
                case 'M':
                    padded_output.push_back(sequence[incoming_index]);
                    if(!quality.empty())
                        padded_quality.push_back(quality[incoming_index]);

                    incoming_index += 1;
                    template_index += 1;
                    cigar_index += 1;
                    break;
                case 'D':
                    insertGapBeforeColumn(template_index + template_leading);
                    padded_output.push_back(sequence[incoming_index]);
                    if(!quality.empty())
                        padded_quality.push_back(quality[incoming_index]);

                    incoming_index += 1;
                    cigar_index += 1;
                    template_index += 1; // skip the newly introduced gap
                    break;
                case 'I':
                    padded_output.push_back('-');
                    if(!quality.empty())
                        padded_quality.push_back('-');
                        
                    cigar_index += 1;
                    template_index += 1;
                    break;
            }
        }
    }

    // Calculate the number of unfilled columns of the multiple alignment that come after
    // the padded sequence
    size_t incoming_trailing = template_element->getNumColumns() - padded_output.size() - incoming_leading;
    MultipleAlignmentElement incoming_element(name, padded_output, padded_quality, 
                                              incoming_leading, incoming_trailing);

    m_sequences.push_back(incoming_element);
}

//
void MultipleAlignment::print(int max_columns) const
{
    (void)max_columns;
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        std::string padding = std::string(m_sequences[i].leading_columns, ' ');
        printf("\t%s%s\t%s\n", padding.c_str(),
                               m_sequences[i].padded_sequence.c_str(), 
                               m_sequences[i].name.c_str());
    }
}

//
void MultipleAlignment::printPileup() const
{
    if(m_sequences.empty())
        return;

    // Get the total number of columns in the alignment
    size_t num_columns = m_sequences.front().getNumColumns();
    size_t num_sequences = m_sequences.size();
    for(size_t i = 0; i < num_columns; ++i) {
        std::string pileup;
        std::string quality;
        for(size_t j = 0; j < num_sequences; ++j) {
            // sanity check that the columns are set up correctly
            assert(m_sequences[j].getNumColumns() == num_columns);
            char symbol = m_sequences[j].getColumnSymbol(i);
            if(symbol != '\0')
                pileup.push_back(symbol);

            char quality_symbol =  m_sequences[j].getColumnQuality(i);
            if(quality_symbol != '\0')
                quality.push_back(quality_symbol);

        }
        printf("%zu\t%s\t%s\n", i, pileup.c_str(), quality.c_str());
    }
}

//
void MultipleAlignment::insertGapBeforeColumn(size_t column_index)
{
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        m_sequences[i].insertGapBeforeColumn(column_index);
    }
}

//
std::string MultipleAlignment::expandCigar(const std::string& cigar)
{
    std::string out;
    std::stringstream parser(cigar);
    int length;
    char symbol;
    while(parser >> length >> symbol)
        out.append(length, symbol);
    return out;
}
