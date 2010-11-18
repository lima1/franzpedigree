/*
 * $Id: io_tests.c 1672 2009-06-18 15:16:21Z markus $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tap.h"
#include "global.h"
#include "dataio.h"
#include "prob.h"
#include "macros.h"
#include "freq.h"

DATA Data;
PROBS Probs;
OPTIONS Options;

int
main(int argc, char *argv[])
{
    int i, ***gcloned;
    char* filename   = FILENAME;
    char* filenameMD = FILENAMEMD;
    char* filenamePOPS = FILENAMEPOPS;
    char* filenamePOPSDIST = FILENAMEPOPSDIST;
    char* filenameMOTHERS = FILENAMEMOTHERS;

    plan_tests(55);

    DATAIOinfile(filename);
    strncpy(Options.PedigreeInfilename,filenameMOTHERS, PATHLEN);
    ok(Data.num_samples == 5, "5 samples in test data");
    gcloned = DATAIOcloneGenotypes();
    ok(gcloned[4][1][1] == 188, "cloning works");
    DATAIOdestroyGenotypes(gcloned);
    DATAIOdestroyPreMCMC();
    DATAIOdestroyPostMCMC();
    Options.MinTyped = 2;
    DATAIOinfile(filenameMD);
    DATAIOinitInpedigree();

    ok(Data.num_samples == 7, "7 samples in simpsons data");
    ok(Data.num_loci == 3, "3 loci in simpsons data");
    ok(Data.locus_has_missing_data[0] == 0, "First locus complete");
    ok(Data.locus_has_missing_data[1] == 0, "Second locus complete");
    ok(Data.locus_has_missing_data[2] == 1, "Third locus incomplete");
    ok(strcmp(Data.samples[0][0].description,"Grampa") == 0, "Description correct");
    ok(strcmp(Data.samples[0][1].description,"Homer") == 0, "Description correct");
    ok(strcmp(Data.samples[0][2].description,"Bart") == 0, "Description correct");
    ok(strcmp(Data.samples[0][3].description,"Lisa") == 0, "Description correct");
    ok(strcmp(Data.samples[0][4].description,"Maggie3210") == 0, "Description correct");
    ok(strcmp(Data.samples[0][5].description,"Marge") == 0, "Description correct");
    ok(strcmp(Data.samples[0][6].description,"Flanders") == 0, "Description correct");
    ok(Data.samples[0][0].typed_loci == 2, "Typed loci correct");
    ok(Data.samples[0][1].typed_loci == 3, "Typed loci correct");
    ok(Data.samples[0][2].typed_loci == 3, "Typed loci correct");
    ok(Data.samples[0][3].typed_loci == 3, "Typed loci correct");
    ok(Data.samples[0][4].typed_loci == 3, "Typed loci correct");
    ok(Data.samples[0][5].typed_loci == 2, "Typed loci correct");
    ok(Data.samples[0][6].typed_loci == 2, "Typed loci correct");
    ok(!HAS_KNOWN_MOTHER(0), "no known mother");
    ok(!HAS_KNOWN_MOTHER(1), "no known mother");
    ok(HAS_KNOWN_MOTHER(2), "known mother");
    ok(HAS_KNOWN_MOTHER(3), "known mother");
    ok(HAS_KNOWN_MOTHER(4), "known mother");
    ok(!HAS_KNOWN_MOTHER(5), "no known mother");
    ok(!HAS_KNOWN_MOTHER(6), "no known mother");

    ok(Data._num_missing_data == 4, "4 alleles untyped");
    ok(Data._missing_data[0].both && Data._missing_data[1].both, "grampa two untyped allele");
    ok(Data._missing_data[0].id == 0  && Data._missing_data[1].id == 0, "grampa two untyped allele");
    ok(Data._missing_data[0].locus == 2  && Data._missing_data[1].locus == 2, "grampa two untyped allele at locus 2");
    ok(Data._missing_data[0].allele == 0  && Data._missing_data[1].allele == 1, "grampa two untyped allele at locus 2");
    ok(Data._missing_data[2].id == 5,  "marge one missing allele");
    ok(Data._missing_data[3].id == 6,  "flanders one missing allele");
    ok(Data._missing_data[2].allele == 0,  "marge one missing allele");
    ok(Data._missing_data[3].allele == 0,  "flanders one missing allele");
    ok(Data._missing_data[2].locus == 2  && Data._missing_data[3].locus == 2, "marge/flanders at locus 2");
    
    
    DATAIOdestroyPreMCMC();
    DATAIOdestroyPostMCMC();
    DATAIOinfile(filenamePOPS);
    ok(Data.num_populations == 3, "3 populations in data");
    ok(!Data.use_distances, "no distances yet");
    DATAIOgeofile(filenamePOPSDIST, 1);
    ok(Data.use_distances, "distances available");

    for (i=0; i<Data.num_populations; i++)
        ok(CMP_DBL(Data.population_distances[i][i], 0.), "zero distance");
    
    ok(CMP_DBL(Data.population_distances[0][1], 0.), "distance correct");
    ok(fabs(Data.population_distances[0][2] - 1030.116) < 0.001, "distance correct");
    ok(fabs(Data.population_distances[2][0] - 1030.116) < 0.001, "distance correct");
    ok(fabs(Data.population_distances[2][1] - 1030.116) < 0.001, "distance correct");
    DATAIOgeofile(FILENAMEPOPSCOORD, 0);
    for (i=0; i<Data.num_populations; i++)
        ok(CMP_DBL(Data.population_distances[i][i],0.), "zero distance");
    
    ok(CMP_DBL(Data.population_distances[0][1], 0.), "distance correct");
    ok(fabs((Data.population_distances[0][2]/1000.) - 1030.116) < 0.001, "distance correct");
    ok(fabs((Data.population_distances[2][0]/1000.) - 1030.116) < 0.001, "distance correct");
    ok(fabs((Data.population_distances[2][1]/1000.) - 1030.116) < 0.001, "distance correct");
    
    DATAIOdestroyPreMCMC();
    DATAIOdestroyPostMCMC();
    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
