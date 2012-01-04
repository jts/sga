//-------------------------------------------------------------------------------
// 
// overlapper - String-string overlap algorithm 
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
#include "overlapper.h"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <limits>
#include <stdio.h>

//
#define max3(x,y,z) std::max(std::max(x,y), z)
//#define DEBUG_OVERLAPPER 1

//
bool SequenceOverlap::isValid() const
{
    return !cigar.empty() && match[0].isValid() && match[1].isValid();
}

//
double SequenceOverlap::getPercentIdentity() const
{
    return (double)(total_columns - edit_distance) * 100.0f / total_columns;
}

//
std::ostream& operator<<(std::ostream& out, const SequenceOverlap& overlap)
{
    out << "[" << overlap.match[0].start << " " << overlap.match[0].end << "] ";
    out << "[" << overlap.match[1].start << " " << overlap.match[1].end << "] ";
    out << "C:" << overlap.cigar;
    return out;
}

//
void SequenceOverlap::printAlignment(const std::string& s1, const std::string& s2) const
{
    assert(isValid());

    std::string out_1;
    std::string out_2;

    // Print out the initial part of the strings, which do not match. 
    // Typically this is the overhanging portion of one of the strings.
    std::string leader_1 = s1.substr(0, match[0].start);
    std::string leader_2 = s2.substr(0, match[1].start);

    // Pad the beginning of the output strings with spaces to align
    if(leader_1.size() < leader_2.size())
        out_1.append(leader_2.size() - leader_1.size(), ' ');

    if(leader_2.size() < leader_1.size())
        out_2.append(leader_1.size() - leader_2.size(), ' ');
    
    out_1.append(leader_1);
    out_2.append(leader_2);

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(s2.substr(current_2, length));
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'I') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(length, '-');
            current_1 += length;
        }
        else if(code == 'D') {
            out_1.append(length, '-');
            out_2.append(s2.substr(current_2, length));
            current_2 += length;
        }
        length = -1;
    }

    // Append the remainder of each string
    out_1.append(s1.substr(current_1));
    out_2.append(s2.substr(current_2));

    // Print the output strings and split long lines
    int MAX_COLUMNS = 120;
    size_t total_columns = std::max(out_1.size(), out_2.size());
    for(size_t i = 0; i < total_columns; i += MAX_COLUMNS) {
        std::string sub_1;
        std::string sub_2;
        if(i < out_1.size())
            sub_1 = out_1.substr(i, MAX_COLUMNS);
        if(i < out_2.size())
            sub_2 = out_2.substr(i, MAX_COLUMNS);
        
        std::cout << "S1\t" << sub_1 << "\n";
        std::cout << "S2\t" << sub_2 << "\n";
        std::cout << "\n";
    }
    std::cout << "Score: " << score << "\n";
    printf("Identity: %2.2lf\n", getPercentIdentity());
}

//
SequenceOverlap Overlapper::computeOverlap(const std::string& s1, const std::string& s2)
{
    // Exit with invalid intervals if either string is zero length
    SequenceOverlap output;
    if(s1.empty() || s2.empty()) {
        std::cerr << "Overlapper::computeOverlap error: empty input sequence\n";
        exit(EXIT_FAILURE);
    }

    // We use same scoring as bwasw
    const int MATCH_SCORE = 2;
    const int GAP_PENALTY = -5;
    const int MISMATCH_PENALTY = -3;

    // Initialize the scoring matrix
    size_t num_columns = s1.size() + 1;
    size_t num_rows = s2.size() + 1;

    typedef std::vector<int> IntVector;
    std::vector<IntVector> score_matrix;
    score_matrix.resize(num_columns);
    for(size_t i = 0; i < score_matrix.size(); ++i)
        score_matrix[i].resize(num_rows);

    // Calculate scores
    for(size_t i = 1; i < num_columns; ++i) {
        for(size_t j = 1; j < num_rows; ++j) {
            // Calculate the score for entry (i,j)
            int idx_1 = i - 1;
            int idx_2 = j - 1;
            int diagonal = score_matrix[i-1][j-1] + (s1[idx_1] == s2[idx_2] ? MATCH_SCORE : MISMATCH_PENALTY);
            int up = score_matrix[i][j-1] + GAP_PENALTY;
            int left = score_matrix[i-1][j] + GAP_PENALTY;

            score_matrix[i][j] = max3(diagonal, up, left);
        }
    }
 
    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(size_t i = 1; i < num_columns; ++i) {
        int v = score_matrix[i][num_rows - 1];
        if(score_matrix[i][num_rows - 1] > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(size_t j = 1; j < num_rows; ++j) {
        int v = score_matrix[num_columns - 1][j];
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }

    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_OVERLAPPER
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];

        int diagonal = score_matrix[i - 1][j - 1] + (is_match ? MATCH_SCORE : MISMATCH_PENALTY);
        int up = score_matrix[i][j-1] + GAP_PENALTY;
        int left = score_matrix[i-1][j] + GAP_PENALTY;

        if(score_matrix[i][j] == diagonal) {
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }
        else if(score_matrix[i][j] == up) {
            cigar.push_back('D');
            j -= 1;
            output.edit_distance += 1;
        }
        else {
            assert(score_matrix[i][j] == left);
            cigar.push_back('I');
            i -= 1;
            output.edit_distance += 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());

    std::stringstream compact_cigar;
    char curr_symbol = cigar[0];
    int curr_run = 1;
    for(size_t i = 1; i < cigar.size(); ++i) {
        if(cigar[i] == curr_symbol) {
            curr_run += 1;
        }
        else {
            compact_cigar << curr_run << curr_symbol;
            curr_symbol = cigar[i];
            curr_run = 1;
        }
    }

    // Add last symbol/run
    compact_cigar << curr_run << curr_symbol;
    output.cigar = compact_cigar.str();
    return output;
}

