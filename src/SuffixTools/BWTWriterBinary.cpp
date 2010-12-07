//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Write a run-length encoded BWT to a binary file
//
#include "BWTWriterBinary.h"
#include "SBWT.h"
#include "RLBWT.h"

//
BWTWriterBinary::BWTWriterBinary(const std::string& filename) : m_numRuns(0), m_runFileOffset(0), m_stage(IOS_NONE)
{
    m_pWriter = createWriter(filename, std::ios::out | std::ios::binary);
    m_stage = IOS_HEADER;
}

//
BWTWriterBinary::~BWTWriterBinary()
{
    assert(m_stage == IOS_DONE);
    delete m_pWriter;
}

//
void BWTWriterBinary::writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag)
{
    size_t largeSampleRate = RLBWT::DEFAULT_SAMPLE_RATE_LARGE;
    size_t smallSampleRate = RLBWT::DEFAULT_SAMPLE_RATE_SMALL;
    m_compressor.initialize(largeSampleRate, smallSampleRate, num_symbols);

    assert(m_stage == IOS_HEADER);
    m_pWriter->write(reinterpret_cast<const char*>(&RLBWT_FILE_MAGIC), sizeof(RLBWT_FILE_MAGIC));
    m_pWriter->write(reinterpret_cast<const char*>(&num_strings), sizeof(num_strings));
    m_pWriter->write(reinterpret_cast<const char*>(&num_symbols), sizeof(num_symbols));

    // Write the sample rates
    m_pWriter->write(reinterpret_cast<const char*>(&largeSampleRate), sizeof(largeSampleRate));
    m_pWriter->write(reinterpret_cast<const char*>(&smallSampleRate), sizeof(smallSampleRate));

    // Here we do not know the number of bytes that are going to be written to the file
    // for the compressed bwt string so we save the offset in the file and write a dummy value. 
    //After the bwt string has been written, we return here and fill in the correct value
    m_runFileOffset = m_pWriter->tellp();

    size_t tmpUnits = 0;
    m_pWriter->write(reinterpret_cast<const char*>(&tmpUnits), sizeof(tmpUnits));

    assert(flag == BWF_NOFMI);
    m_pWriter->write(reinterpret_cast<const char*>(&flag), sizeof(flag));

    m_stage = IOS_BWSTR;    
}

// Write a single character of the BWStr
// If the char is '\n' we are finished
void BWTWriterBinary::writeBWChar(char b)
{
    m_compressor.writeSymbol(b, m_pWriter);
}

// write the final run to the stream and fill in the number of runs
void BWTWriterBinary::finalize()
{
    m_compressor.flush(m_pWriter);

    size_t numStringBytes = m_compressor.getNumBytesWrote();
    size_t savedPos = m_pWriter->tellp();
    m_pWriter->seekp(m_runFileOffset);
    m_pWriter->write(reinterpret_cast<const char*>(&numStringBytes), sizeof(numStringBytes));
    
    // Return to the end of the file
    m_pWriter->seekp(savedPos);

    // Write out the large and small markers
    const LargeMarkerVector& largeMarkers = m_compressor.getLargeMarkerVector();
    size_t numLargeMarkers = largeMarkers.size();

    m_pWriter->write(reinterpret_cast<const char*>(&numLargeMarkers), sizeof(numLargeMarkers));
    for(size_t i =  0; i < numLargeMarkers; ++i)
        m_pWriter->write(reinterpret_cast<const char*>(&largeMarkers[i]), sizeof(LargeMarker));

    const SmallMarkerVector& smallMarkers = m_compressor.getSmallMarkerVector();
    size_t numSmallMarkers = smallMarkers.size();

    m_pWriter->write(reinterpret_cast<const char*>(&numSmallMarkers), sizeof(numSmallMarkers));
    for(size_t i =  0; i < numSmallMarkers; ++i)
        m_pWriter->write(reinterpret_cast<const char*>(&smallMarkers[i]), sizeof(SmallMarker));

    // Write out the Pred array
    AlphaCount64 finalCounts = m_compressor.getRunningCount();
    
    AlphaCount64 predCount;
    predCount.set('$', 0);
    predCount.set('A', finalCounts.get('$')); 
    predCount.set('C', predCount.get('A') + finalCounts.get('A'));
    predCount.set('G', predCount.get('C') + finalCounts.get('C'));
    predCount.set('T', predCount.get('G') + finalCounts.get('G'));
    m_pWriter->write(reinterpret_cast<const char*>(&predCount), sizeof(predCount));

    m_stage = IOS_DONE;
}

