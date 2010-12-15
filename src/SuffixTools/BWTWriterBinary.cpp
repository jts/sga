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
BWTWriterBinary::BWTWriterBinary(const std::string& filename, 
                                 int smallSampleRate) : m_smallSampleRate(smallSampleRate),
                                                          m_filename(filename), 
                                                          m_numRuns(0), 
                                                          m_headerFileOffset(0), 
                                                          m_stage(IOS_NONE)
{
    m_pWriter = createWriter(m_filename, std::ios::out | std::ios::binary);
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
    m_compressor.initialize(m_filename, largeSampleRate, m_smallSampleRate, num_symbols);

    assert(m_stage == IOS_HEADER);
    m_pWriter->write(reinterpret_cast<const char*>(&RLBWT_FILE_MAGIC), sizeof(RLBWT_FILE_MAGIC));
    m_pWriter->write(reinterpret_cast<const char*>(&num_strings), sizeof(num_strings));
    m_pWriter->write(reinterpret_cast<const char*>(&num_symbols), sizeof(num_symbols));

    // Write the sample rates
    m_pWriter->write(reinterpret_cast<const char*>(&largeSampleRate), sizeof(largeSampleRate));
    m_pWriter->write(reinterpret_cast<const char*>(&m_smallSampleRate), sizeof(m_smallSampleRate));

    // Here we do not know the number of units or the number of markers that are going to be written to the file.
    // We save the offset of this position in the file so we can return later to write in these values.
    m_headerFileOffset = m_pWriter->tellp();

    // Write placeholders for the number of units in the bwt string
    // and the position in the file of the small and large marker segments
    size_t tmpUnits = 0;
    m_pWriter->write(reinterpret_cast<const char*>(&tmpUnits), sizeof(tmpUnits));
    m_pWriter->write(reinterpret_cast<const char*>(&tmpUnits), sizeof(tmpUnits));
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

    // Write out the large and small markers
    size_t largeMarkersStart = m_pWriter->tellp();
    m_compressor.writeLargeMarkers(m_pWriter);

    size_t smallMarkersStart = m_pWriter->tellp();
    m_compressor.writeSmallMarkers(m_pWriter);

    // Write out the Pred array
    AlphaCount64 finalCounts = m_compressor.getRunningCount();
    
    AlphaCount64 predCount;
    predCount.set('$', 0);
    predCount.set('A', finalCounts.get('$')); 
    predCount.set('C', predCount.get('A') + finalCounts.get('A'));
    predCount.set('G', predCount.get('C') + finalCounts.get('C'));
    predCount.set('T', predCount.get('G') + finalCounts.get('G'));
    m_pWriter->write(reinterpret_cast<const char*>(&predCount), sizeof(predCount));

    // Return to the header of the file and overwrite the placeholders with the real values
    size_t numStringBytes = m_compressor.getNumBytesWrote();
    size_t savedPos = m_pWriter->tellp();
    m_pWriter->seekp(m_headerFileOffset);
    m_pWriter->write(reinterpret_cast<const char*>(&numStringBytes), sizeof(numStringBytes));
    m_pWriter->write(reinterpret_cast<const char*>(&largeMarkersStart), sizeof(largeMarkersStart));
    m_pWriter->write(reinterpret_cast<const char*>(&smallMarkersStart), sizeof(smallMarkersStart));
    
    // Return to the end of the file
    m_pWriter->seekp(savedPos);

    m_stage = IOS_DONE;
}

