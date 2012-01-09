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
//#define MA_DEBUG_CONSENSUS 1

// Initialize static members
const char* MultipleAlignment::m_alphabet = "ACGTN-";

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
size_t MultipleAlignmentElement::getStartColumn() const
{
    return leading_columns;
}

//
size_t MultipleAlignmentElement::getEndColumn() const
{
    return getNumColumns() - trailing_columns - 1;
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
std::string MultipleAlignmentElement::getPrintableSubstring(size_t start_column, size_t num_columns) const
{
    std::string out;
    size_t i = 0;
    while(i < num_columns) {
        char symbol = getColumnSymbol(start_column + i);
        out.push_back(symbol != '\0' ? symbol : ' ');
        ++i;
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
            if(expanded_cigar[cigar_index] == 'I') {
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
                case 'I':
                    insertGapBeforeColumn(template_index + template_leading);
                    padded_output.push_back(sequence[incoming_index]);
                    if(!quality.empty())
                        padded_quality.push_back(quality[incoming_index]);

                    incoming_index += 1;
                    cigar_index += 1;
                    template_index += 1; // skip the newly introduced gap
                    break;
                case 'D':
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

std::string MultipleAlignment::calculateBaseConsensus(int min_call_coverage, int min_trim_coverage)
{
    assert(!m_sequences.empty());
    std::string consensus_sequence;
    MultipleAlignmentElement& base_element = m_sequences.front();
    size_t start_column = base_element.getStartColumn();
    size_t end_column = base_element.getEndColumn();
    
    // This index records the last base in the consensus that had coverage greater than
    // min_trim_coverage. After the consensus calculation the read is trimmed back to this position
    int last_good_base = -1;

    for(size_t c = start_column; c <= end_column; ++c) {
        std::vector<int> counts = getColumnBaseCounts(c);

        char max_symbol = '\0';
        int max_count = -1;
        int total_depth = 0;

#ifdef MA_DEBUG_CONSENSUS
        printf("%zu\t", c);
#endif
        for(size_t a = 0; a < m_alphabet_size; ++a) {
            char symbol = m_alphabet[a];
            total_depth += counts[a];

            if(symbol != 'N' && counts[a] > max_count) {
                max_symbol = symbol;
                max_count = counts[a];
            }
#ifdef MA_DEBUG_CONSENSUS
            printf("%c:%d ", symbol, counts[a]);
#endif
        }

        char base_symbol = base_element.getColumnSymbol(c);
        int base_count = counts[symbol2index(base_symbol)];
        
        // Choose a consensus base for this column. Only change a base
        // if has been seen less than min_call_coverage times and the max
        // base in the column has been seen more times than the base symbol
        char consensus_symbol;
        if(max_count > base_count && base_count < min_call_coverage)
            consensus_symbol = max_symbol;
        else
            consensus_symbol = base_symbol;
        
        // Output a symbol to the consensus. Skip padding symbols and leading
        // bases that are less than the minimum required depth to avoid trimming
        if(consensus_symbol != '-' &&
            (!consensus_sequence.empty() || total_depth >= min_trim_coverage))
                consensus_sequence.push_back(consensus_symbol);

        // Record the position of the last good base
        if(total_depth >= min_trim_coverage) {
            int consensus_index = consensus_sequence.size() - 1;
            if(consensus_index > last_good_base)
                last_good_base = consensus_index;
        }
#ifdef MA_DEBUG_CONSENSUS
        printf("CALL: %c\n", consensus_symbol);
#endif 
    }

    if(last_good_base != -1)
        consensus_sequence.erase(last_good_base);
    else
        consensus_sequence.clear();

    return consensus_sequence;
}

//
void MultipleAlignment::print(size_t max_columns) const
{
    if(m_sequences.empty())
        return;

    size_t total_columns = m_sequences.front().getNumColumns();

    // Print the multiple alignment in segments
    for(size_t c = 0; c < total_columns; c += max_columns) {
        size_t remaining = total_columns - c;
        size_t slice_size = max_columns < remaining ? max_columns : remaining;
        for(size_t i = 0; i < m_sequences.size(); ++i) {
            std::string slice =  m_sequences[i].getPrintableSubstring(c, slice_size);

            // Check if this string is blank, if so don't print it
            if(slice.find_first_not_of(" ") != std::string::npos)
                printf("\t%s\t%s\n", slice.c_str(), m_sequences[i].name.c_str());
        }
        printf("\n\n");
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
        std::string counts_str = getColumnCountString(i);
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
        printf("%zu\t%s\t%s\t%s\n", i, pileup.c_str(), quality.c_str(), counts_str.c_str());
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

//
int MultipleAlignment::symbol2index(char symbol) const
{
    switch(symbol) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return 4;
        case '-':
            return 5;
    }

    std::cerr << "Error: Unrecognized symbol in multiple alignment\n";
    exit(EXIT_FAILURE);
    return -1;
}

//
std::vector<int> MultipleAlignment::getColumnBaseCounts(size_t idx) const
{
    std::vector<int> out(m_alphabet_size, 0);
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        char symbol = m_sequences[i].getColumnSymbol(idx);
        if(symbol != '\0')
            out[symbol2index(symbol)] += 1;
    }
    return out;
}

std::string MultipleAlignment::getColumnCountString(size_t column) const
{
    std::vector<int> counts = getColumnBaseCounts(column);
    std::stringstream out;
    for(size_t i = 0; i < m_alphabet_size; ++i) {
        out << m_alphabet[i] << ":" << counts[i] << " ";
    }
    return out.str();
}


