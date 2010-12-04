//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#include <iostream>
#include <math.h>
#include "Util.h"

//
// Sequence operations
//

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t last_pos = seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
    {
        out[last_pos - i] = complement(seq[i]);
    }
    return out;
}

// Reverse a sequence
std::string reverse(const std::string& seq)
{
    return std::string(seq.rbegin(), seq.rend());
}

// Complement a sequence
std::string complement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t l = seq.length();
    for(size_t i = 0; i < l; ++i)
    {
        out[i] = complement(seq[i]);
    }
    return out;
}

// 
std::string prefix(const std::string& seq, const unsigned int len)
{
    assert(seq.length() >= len);
    return seq.substr(0, len);
}

// 
std::string suffix(const std::string& seq, const unsigned int len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}

// count the differences between s1 and s2 over the first n bases
int countDifferences(const std::string& s1, const std::string& s2, size_t n)
{
    int numDiff = 0;
    for(size_t i = 0; i < n; ++i)
    {
        if(s1[i] != s2[i])
            numDiff++;
    }
    return numDiff;
}

// count the differences between s1 and s2 over the first n bases
std::string getDiffString(const std::string& s1, const std::string& s2)
{
    std::string out;
    size_t stop = std::min(s1.size(), s2.size());
    for(size_t i = 0; i < stop; ++i)
        out.push_back(s1[i] == s2[i] ? '.' : s1[i]);
    return out;
}

// 
char randomBase()
{
    int i = rand() % 4;
    switch(i)
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            assert(false);
    }
    return 'A';        
}

// Strip the leading directories and
// the last trailling suffix from a filename
std::string stripFilename(const std::string& filename)
{
    std::string out = stripDirectories(filename);
    // Remove the gzip extension if necessary
    if(isGzip(out))
        out = stripExtension(out);
    return stripExtension(out);
}

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}

// Strip the leadering directories from a filename
std::string stripDirectories(const std::string& filename)
{
    size_t lastDirPos = filename.find_last_of('/');
    
    if(lastDirPos == std::string::npos)
        return filename; // no directories
    else
        return filename.substr(lastDirPos + 1);
}

// Return the file extension
std::string getFileExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    
    if(suffixPos == std::string::npos)
        return "";
    else
        return filename.substr(suffixPos + 1);
}

// slow function to convert an integer to a binary string
std::string int2Binary(size_t v, int numBits)
{
    std::string tmp;
    int bits = sizeof(v) * 8;
    for(int i = bits - 1; i >= 0; --i)
    {
        // test if the i-th bit is set
        size_t mask = 1;
        mask <<= i;
        char b = (v & mask) ? '1' : '0';
        tmp.append(1, b);
    }

    size_t pos = 0;
    if(numBits == 0)
    {
        // Truncate leading zeros
        size_t pos = tmp.find_first_of('1');
        if(pos == std::string::npos)
            pos = tmp.size() - 1;
    }
    else
    {
        pos = tmp.size() - numBits;
    }
    return tmp.substr(pos);
}

// Write out a fasta record
void writeFastaRecord(std::ostream* pWriter, const std::string& id, const std::string& seq, size_t maxLineLength)
{
    *pWriter << ">" << id << " " << seq.length() << "\n";
    for(size_t i = 0; i < seq.size(); i += maxLineLength)
    {
        size_t span = std::min(seq.size() - i, maxLineLength);
        *pWriter << seq.substr(i, span) << "\n";
    }
}

// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    std::string extension = suffix(filename, sizeof(GZIP_EXT) - 1);
    return extension == GZIP_EXT;
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        ogzstream* pGZ = new ogzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
        assertFileOpen(*pWriter, filename);
        return pWriter;
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }    
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }    
}

//
void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Split a string into parts based on the delimiter
StringVector split(std::string in, char delimiter)
{
    StringVector out;
    size_t lastPos = 0;
    size_t pos = in.find_first_of(delimiter);

    while(pos != std::string::npos)
    {
        out.push_back(in.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
        pos = in.find_first_of(delimiter, lastPos);
    }
    out.push_back(in.substr(lastPos));
    return out;
}

// Split a key-value pair
void splitKeyValue(std::string in, std::string& key, std::string& value)
{
    StringVector parts = split(in, CAF_SEP);
    if(parts.size() != 2)
    {
        std::cerr << "Invalid key-value pair " << in << std::endl;
        assert(false);
    }

    key = parts[0];
    value = parts[1];

    assert(key.size() > 0 && value.size() > 0 && "Invalid key-value pair");
}

// Get the shared component of a read pair name
// This is the part of the name preceding the "/"
std::string getPairBasename(const std::string& id)
{
    assert(!id.empty());

    size_t pos = id.find_last_of('/');
    if(pos == std::string::npos)
        return id;
    else
        return id.substr(0, pos);
}

// Get the ID of the pair of a given read
std::string getPairID(const std::string& id)
{
    assert(!id.empty());
    std::string pid(id);

    size_t li = id.length() - 1;
    char last = id[li];

    if(last == 'A')
        pid[li] = 'B';
    else if(last == 'B')
        pid[li] = 'A';
    else if(last == '1')
        pid[li] = '2';
    else if(last == '2')
        pid[li] = '1';
    else
        pid = "";
    return pid;
}

// Return 0 if the id indicates the first read of a pair, 1 otherwise
int getPairIndex(const std::string& id)
{
    size_t li = id.length() - 1;
    char last = id[li];
    if(last == 'A' || last == '1')
        return 0;
    else if(last == 'B' || last == '2')
        return 1;
    else
    {
        std::cerr << "Unrecognized pair format: " << id << "\n";
        exit(EXIT_FAILURE);
    }
}

// Debug function to parse the distance between two reads
// based on their names
// This assumes is that the read names are just the positions the reads
// were sampled from
size_t debug_getReadDistFromNames(const std::string& name1, const std::string& name2)
{
    std::string id;
    std::string pos;
    StringVector parts = split(name1, ':');
    if(parts.size() != 2)
        return 0;

    int p1 = atoi(parts[1].c_str());
    
    parts = split(name2, ':');
    if(parts.size() != 2)
        return 0;

    int p2 = atoi(parts[1].c_str());
    int dist = int(abs(p1 - p2));
    return dist;
}

