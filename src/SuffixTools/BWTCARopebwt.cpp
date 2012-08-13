//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTCARopebwt - Construct the BWT for a set of reads
// using Heng Li's ropebwt implementation
//
#include "BWTCARopebwt.h"
#include "bcr.h"
#include "SeqReader.h"
#include "StdAlnTools.h"

void BWTCA::runRopebwt(const std::string& input_filename)
{
    // Initialize ropebwt
    bcr_t* bcr = bcr_init(true, "ropebwt-tmp");

    SeqReader reader(input_filename);
    SeqRecord record;
    while(reader.get(record))
    {
        // Convert the sequence of this read to a 2-bit packed representation
        // ropebwt uses Heng Li's standard 2-bit sequence encoding, which is wrapped in StdalnTools
        uint8_t* s = StdAlnTools::createPacked(record.seq.toString());
        bcr_append(bcr, record.seq.length(), s);
        delete [] s;
    }

    std::cout << "Done reading sequences\n";
    // Build the BWT
    bcr_build(bcr);
    
    std::cout << "Done build\n";

#define print_bwt(itr_t, itr_set, itr_next_f, is_bin, fp) do { \
        itr_t *itr; \
        const uint8_t *s; \
        int i, j, l; \
        itr = (itr_set); \
        if (is_bin) { \
            fwrite("RLE\6", 4, 1, fp); \
            while ((s = itr_next_f(itr, &l)) != 0) \
                fwrite(s, 1, l, fp); \
        } else { \
        while ((s = itr_next_f(itr, &l)) != 0) \
            for (i = 0; i < l; ++i) \
                for (j = 0; j < s[i]>>3; ++j) \
                    fputc("$ACGTN"[s[i]&7], fp); \
                    fputc('\n', fp); \
                } \
            free(itr); \
        } while (0)


    FILE* out = fopen("ropebwt.bwt", "w");
    print_bwt(bcritr_t, bcr_itr_init(bcr), bcr_itr_next, 1, out);
    bcr_destroy(bcr);
}
