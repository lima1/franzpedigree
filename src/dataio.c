/*
 * $Id: dataio.c 2062 2010-05-25 10:24:32Z markus $
 *
 * Copyright (C) 2008-2010. Universitaet Leipzig  
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.*
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include "macros.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <ctype.h>              /* for toupper */

#include "dataio.h"

#include "options.h"
#include "prob.h"
#include "freq.h"
#include "global.h"
#include "dag.h"
#include "uthash.h"
#include "utils.h"
#include "lod.h"

extern DATA Data;
extern PROBS Probs;
extern OPTIONS Options;

/* read the next line of an input file                                 */
#define FGETS if ((fgets(input, (int)sizeof(input), in)) == NULL)            \
    fatalParseError(filename, msg_prefix, "truncated inputfile", line); \
    line++;                                                             \

static void
fatalParseError(char *filename, char *msg_prefix, char *msg, int line)
{
    int help_id;
    char message[LINESIZE];
    char help[][LINESIZE] = { "1 3 / SIMPSONS\n^ ^ ^ ^\n| | | |\n| | | Description\n| | Separator\n| Number Loci\nNumber Populations/Locations\n",
                              "7 Springfield\n^ ^\n| |\n| Description\nNumber genotypes\n\n (see manual if you want to provide loci ids)\n",
                              "Grampa     1 1920 ? M 110/100 200/208 ?/?\n|          | |    | | |\n|          | |    | | Alleles\n|          | |    | Sex (M, F or ?)\n|          | |    Year of death (? = unkown)\n|          | Year of birth\n|         Ramets\nDescription (exactly 10 characters)\n\nMOST COMMON ERROR: Description too short or too long. Must be exactly ten characters long!\n"
 };
    switch (line) {
        case 1:
            help_id = 0;
            break;
        case 2:
            help_id = 1;
            break;
        default:
            help_id = 2;
    }
    if (strncmp(msg_prefix, "INFILE", 6) == 0 )
        snprintf(message, sizeof(message),
                "[while parsing file %s,line %i] %s%s\n\nExample of a valid line %i:\n\n%s", filename, line,
                msg_prefix, msg, line, help[help_id]);
    else 
        snprintf(message, sizeof(message),
                "[while parsing file %s,line %i] %s%s\n", filename, line, msg_prefix, msg);

    FATAL(message);
}

void
DATAIOinitInpedigree()
{
    int i, j;

    MAKE2DINT(Data.in_pedigree, Data.num_samples, 2, i);
    MAKE1DINT(Data.in_pedigree_relationships, Data.num_loci);
    /* we can use the known relationships to estimate errors   */
    MAKE3DINT(Data.in_pedigree_mismatches, Data.num_samples, 2, Data.num_loci,
              i, j);
    DATAIOpedigreein(Options.PedigreeInfilename, Data.in_pedigree);
    Data.use_pedigree = true;
}

static int
addSampleIDtoHash(int id, char *name)
{
    DESCRIPTION_HASH *s;

    /* we do not yet know the upper bound of id */
    if (id < 0)
        FATALINT("addSampleIDtoHash()");

    HASH_FIND_STR(Data.sample_ids, name, s);
    if (s)
        return 0;

    s = malloc(sizeof(DESCRIPTION_HASH));
    if (s == NULL)
        FATAL("malloc failed");

    s->id = id;
    strcpy(s->name, name);
    HASH_ADD_STR(Data.sample_ids, name, s);
    return 1;
}

static int
getSampleIDfromHash(char *name, char *filename, char *msg_prefix, int line)
{
    DESCRIPTION_HASH *s;
    char msg[LINESIZE];

    HASH_FIND_STR(Data.sample_ids, name, s);
    if (!s) {
        snprintf(msg, sizeof(msg), "unknown id \"%s\"", name);
        fatalParseError(filename, msg_prefix, msg, line);
    }
    /* make sure the id is in range */
    if (s->id < 0 || s->id >= Data.num_samples)
        FATALINT("getSampleIDfromHash()");

    return s->id;
}

void
DATAIOpedigreein(char *filename, int **in_pedigree)
{
    int i, num_samples, line = 0, res, id1, id2;
    FILE *in = NULL;
    char input[LINESIZE], desc1[D_PARSE_DESC_LENGTH+1],
        desc2[D_PARSE_DESC_LENGTH+1];
    char msg_prefix[] = "PEDIGREEINFILE; ";

    /* open file */
    if ((in = fopen(filename, "r")) == NULL)
        FATAL("Cannot open pedigree infile!");

    FGETS;

    res = sscanf(input, "%d", &num_samples);
    if (res != 1)
        fatalParseError(filename, msg_prefix, "invalid header.", line);

    if (Data.num_samples != num_samples)
        fatalParseError(filename, msg_prefix,
                        "different number of samples than input", line);

    for (i = 0; i < Data.num_samples; i++) {
        line++;
        if ((fgets(input, (int)sizeof(input), in)) == NULL)
            fatalParseError(filename, msg_prefix, "missing samples", line);

        UTILScopyNonBlanks(input, desc1, 0, D_PARSE_DESC_LENGTH);
        /* make sure we know this individual */
        (void)getSampleIDfromHash(desc1, filename, msg_prefix, line);
    }
    for (i = 0; i < Data.num_samples; i++) {
        in_pedigree[i][0] = -1;
        in_pedigree[i][1] = -1;
    }
    /* now read the arcs */
    while (fgets(input, (int)sizeof(input), in) != NULL) {
        line++;
        UTILScopyNonBlanks(input, desc1, 0, D_PARSE_DESC_LENGTH);
        id1 = getSampleIDfromHash(desc1, filename, msg_prefix, line);
        UTILScopyNonBlanks(input, desc2, D_PARSE_DESC_LENGTH,
                           D_PARSE_DESC_LENGTH * 2);
        id2 = getSampleIDfromHash(desc2, filename, msg_prefix, line);
        /* id2 already in pedigree? then use second matrix element */
        if (in_pedigree[id2][0] < 0)
            in_pedigree[id2][0] = id1;
        else
            in_pedigree[id2][1] = id1;
    }
    FCLOSE(in);
}

void
DATAIOfullsibFile(char *filename)
{
    int i, j, k, line = 0, res, groups, siblings, *sibling_ids;

    FILE *in = NULL;
    char input[LINESIZE], desc[D_PARSE_DESC_LENGTH+1];
    char msg_prefix[] = "FULLSIBFILE; ";

    if ((in = fopen(filename, "r")) == NULL)
        FATAL("Cannot open fullsib file!");

    DATAIOinitFullsibMatrix();

    FGETS;
    res = sscanf(input, "%d", &groups);
    if (res != 1)
        fatalParseError(filename, msg_prefix,
                        "header lists not the number of fullsib groups",
                        line);

    if (groups < 1 || groups > Data.num_samples / 2)
        fatalParseError(filename, msg_prefix,
                        "number of fullsib groups not in allowed range",
                        line);

    for (i = 0; i < groups; i++) {
        FGETS;
        res = sscanf(input, "%d", &siblings);
        if (res != 1)
            fatalParseError(filename, msg_prefix,
                            "number of siblings not found", line);

        if (siblings < 2 || siblings > Data.num_samples)
            fatalParseError(filename, msg_prefix,
                            "number of siblings not in allowed range", line);
        MAKE1DINT(sibling_ids, siblings);
        for (j = 0; j < siblings; j++) {
            FGETS;
            UTILScopyNonBlanks(input, desc, 0, D_PARSE_DESC_LENGTH);
            sibling_ids[j] =
                getSampleIDfromHash(desc, filename, msg_prefix, line);
        }
        for (j = 0; j < siblings; j++)
            for (k = 0; k < siblings; k++)
                if (Options.FullsibInfile)
                    Data.fs_matrix[TIDXL(sibling_ids[j],sibling_ids[k])].fullsibs =
                        true;
                else  
                    Data.fs_matrix[TIDXL(sibling_ids[j],sibling_ids[k])].halfsibs =
                        true;
        FREE1D(sibling_ids);
    }
}

void
DATAIOfullsibOut(char *filename)
{
    int i,j,cnt_groups = 0,*cnt_ind;
    bool *processed;
    FILE *fp = NULL;

    if (Data.fs_matrix == NULL) FATALINT("DATAIOfullsibOut");
    
    FOPENW(fp, filename);

    /* first we scan one time over the matrix to get the numbers */
    MAKE1DINT(cnt_ind, Data.num_samples);
    MAKE1DBOOL(processed, Data.num_samples);
    for (i=0; i<Data.num_samples; i++) { 
        cnt_ind[i] = 0;
        processed[i] = false;
    }    
    for (i = 0; i < Data.num_samples; i++) {
        if (processed[i]) continue;
        for (j = i+1; j < Data.num_samples; j++)
            if (Data.fs_matrix[TIDX(j,i)].fullsibs) { 
                cnt_ind[i]++;
                processed[j] = true;
            }    
    }
    
    FREE1D(processed);

    for (i=0; i<Data.num_samples; i++) if (cnt_ind[i] > 0) cnt_groups++;

    /* now print the file */
    fprintf(fp, "%d\n", cnt_groups);
    for (i = 0; i < Data.num_samples; i++) {
        if (cnt_ind[i] == 0) continue;
        fprintf(fp, "%d\n%*s\n", (cnt_ind[i]+1), D_PARSE_DESC_LENGTH,Data.id_mapping[i].description);
        for (j = i+1; j < Data.num_samples; j++)
            if (Data.fs_matrix[TIDX(j,i)].fullsibs) 
                fprintf(fp, "%*s\n", D_PARSE_DESC_LENGTH,Data.id_mapping[j].description);
    }
    FREE1D(cnt_ind);
    FCLOSE(fp);
}

void
DATAIOgeofile(char *filename, bool distances)
{
    int i, j, num_populations, line = 0, res;
    FILE *in = NULL;
    char input[LINESIZE], distdata[LINESIZE];
    char msg_prefix[] = "GEO/DISTFILE; ", *p;

    if ((in = fopen(filename, "r")) == NULL)
        FATAL("Cannot open geofile!");

    FGETS;

    Data.max_distance = 0;

    res = sscanf(input, "%d", &num_populations);
    if (res != 1 || Data.num_populations != num_populations)
        fatalParseError(filename, msg_prefix,
                        "different number of populations than input", 1);
    /* already malloc'ed? */
    if (!Data.use_distances)
        MAKE2DDOUBLE(Data.population_distances, num_populations,
                    num_populations, i);
    Data.use_distances = true;

    for (i = 0; i < Data.num_populations; i++) {
        if ((fgets(input, (int)sizeof(input), in) == NULL))
            fatalParseError(filename, msg_prefix,
                            "number of populations probably wrong.", line);
        line++;
        strcpy(distdata, input + D_PARSE_DESC_LENGTH);
        if (distances) {
            p = strtok(distdata, " ");
            if (p == NULL)
                fatalParseError(filename, msg_prefix,
                                "distances wrong (no space separator).", line);
            for (j = 0; j < Data.num_populations; j++) {
                res = sscanf(p, "%lf", &Data.population_distances[i][j]);
                if (res != 1)
                    fatalParseError(filename, msg_prefix,
                                    "could not parse distance", line);
                if (Data.population_distances[i][j] > Data.max_distance)
                    Data.max_distance = Data.population_distances[i][j];

                p = strtok(NULL, " ");
            }
            Data.populations[i].latitude = 0.;
            Data.populations[i].longitude = 0.;
        } else {
            res =
                sscanf(distdata, "%lf %lf", &Data.populations[i].latitude,
                       &Data.populations[i].longitude);
            if (res != 2)
                fatalParseError(filename, msg_prefix,
                                "data not in format latitude longitude",
                                line);
        }
    }
    /* only calc distances if coordinates are now available and user
     * did not provide own distances */
    if (!distances && strlen(Options.Geofilename) == 0) {
        for (i = 0; i < Data.num_populations; i++) {
            Data.population_distances[i][i] = 0;
            for (j = 0; j < Data.num_populations; j++) {
                Data.population_distances[i][j] =
                    UTILScalcDistance(Data.populations[i].latitude,
                                      Data.populations[i].longitude,
                                      Data.populations[j].latitude,
                                      Data.populations[j].longitude);
                if (Data.population_distances[i][j] > Data.max_distance)
                    Data.max_distance = Data.population_distances[i][j];
                Data.population_distances[j][i] =
                    Data.population_distances[i][j];
            }
        }
    }
    FCLOSE(in);
}

void
DATAIOpedigreeout(char *filename, Dag D)
{
    int i, j, indegree;
    FILE *pout;
    Edge incoming[2];

    FOPENW(pout, filename);
    fprintf(pout, "%i\n", Data.num_samples);
    for (i = 0; i < Data.num_samples; i++) {
        fprintf(pout, "%*s\n", D_PARSE_DESC_LENGTH,
                Data.id_mapping[i].description);
    }
    for (i = 0; i < Data.num_samples; i++) {
        indegree = DAGincomingE(D, i, incoming);
        for (j = 0; j < indegree; j++)
            fprintf(pout, "%*s%*s\n", D_PARSE_DESC_LENGTH,
                    Data.id_mapping[incoming[j].v].description,
                    D_PARSE_DESC_LENGTH,
                    Data.id_mapping[incoming[j].w].description);

    }
    FCLOSE(pout);
}

static inline void
cleanID(char *dest, char *src) {
    UTILScopyNonBlanks(src,dest,0, D_PARSE_DESC_LENGTH);
    UTILSreplaceSpace(dest,'_',sizeof(dest));
}

/* output the pedigree in this text format (for the pedigree viewer...) */
void
DATAIOpedigreeoutText(char *filename, Dag D)
{
    int i, indegree;
    FILE *pout;
    char id[2][D_PARSE_DESC_LENGTH+1];
    Edge incoming[2];

    FOPENW(pout, filename);
    fprintf(pout, "%*s %*s %*s\n",
        D_PARSE_DESC_LENGTH, "ID", D_PARSE_DESC_LENGTH, "SIRE", D_PARSE_DESC_LENGTH, "DAM" );
    for (i = 0; i < Data.num_samples; i++) {
        indegree = DAGincomingE(D, i, incoming);
        fprintf(pout, "%*s ", D_PARSE_DESC_LENGTH,Data.id_mapping[i].description);
        if (indegree == 0)
            fprintf(pout, " %*s %*s\n", D_PARSE_DESC_LENGTH, "*", D_PARSE_DESC_LENGTH, "*");
        else if (indegree == 1) {
            cleanID(id[0],Data.id_mapping[incoming[0].v].description);
            if (Data.id_mapping[incoming[0].v].sex == 1) 
                fprintf(pout, " %*s %*s\n", D_PARSE_DESC_LENGTH, "*", D_PARSE_DESC_LENGTH, id[0]);
            else 
                fprintf(pout, " %*s %*s\n", D_PARSE_DESC_LENGTH, id[0], D_PARSE_DESC_LENGTH, "*");
        } else if (indegree == 2) {
            cleanID(id[0],Data.id_mapping[incoming[0].v].description);
            cleanID(id[1],Data.id_mapping[incoming[1].v].description);
            if (Data.id_mapping[incoming[0].v].sex == 1) 
                fprintf(pout, " %*s %*s\n",
                    D_PARSE_DESC_LENGTH, id[1], D_PARSE_DESC_LENGTH, id[0] );
            else 
                fprintf(pout, " %*s %*s\n",
                    D_PARSE_DESC_LENGTH, id[0], D_PARSE_DESC_LENGTH, id[1] );
        }
    }
    FCLOSE(pout);
}

void
DATAIOfreqfile(char *filename, int simulator)
{
    int i, j, num_loci, num_alleles, min, max, allele, line = 0, res;
    double f, fsum;
    FILE *in = NULL;
    char input[LINESIZE], msg_prefix[] = "FREQFILE; ";
    RANGE r;

    /* open file */
    if ((in = fopen(filename, "r")) == NULL)
        FATAL("Cannot open allele frequency file!");

    FGETS;

    res = sscanf(input, "%d", &num_loci);
    /* the simulator uses a freqfile for pedigree simulations */
    if (simulator)
        Data.num_loci = num_loci;
    
    /* in case of --numloci option */
    num_loci = MIN(num_loci, Data.num_loci);

    if (res != 1 || Data.num_loci != num_loci)
        fatalParseError(filename, msg_prefix,
                        "different number of loci than input", line);

    MALLOC(Data.allele_frequencies, FREQS, Data.num_loci);

    for (i = 0; i < num_loci; i++) {
        FGETS;
        res = sscanf(input, "%d %d %d", &num_alleles, &min, &max);
        if (res != 3)
            fatalParseError(filename, msg_prefix,
                            "not in format <num_alleles> <min> <max>", line);

        if (num_alleles < 1 || min < 0 || max < 0 || min > max)
            fatalParseError(filename, msg_prefix,
                            "invalid values for <num_alleles> <min> <max>",
                            line);

        if (num_alleles > D_PARSE_MAX_NUM_ALLELES
            || max - min + 1 > D_PARSE_MAX_ALLELE_RANGE)
            fatalParseError(filename, msg_prefix,
                            "values too large. Change global.h if necessary.",
                            line);

        FREQgetAlleleRange(i, &r);
        MALLOC(Data.allele_frequencies[i].freqs, double, max - min + 1);

        MAKE1DINT(Data.allele_frequencies[i].allele_id_lookup, max - min + 1);

        if (!simulator && (r.max > max || r.min < min))
            fatalParseError(filename, msg_prefix, "allele range wrong", line);
        Data.allele_frequencies[i].min = min;
        Data.allele_frequencies[i].max = max;

        for (j = 0; j < num_alleles; j++) {
            FGETS;
            res = sscanf(input, "%d %lf", &allele, &f);
            if (res != 2)
                fatalParseError(filename, msg_prefix,
                                "not in format <allele> <frequency>", line);
            if (allele < min || allele > max)
                fatalParseError(filename, msg_prefix, "allele not in range.",
                                line);
            if (f < 0. || f > 1.)
                fatalParseError(filename, msg_prefix,
                                "allele frequency not in range.", line);


            Data.allele_frequencies[i].freqs[allele - min] = f;
        }
        fsum = 0.;
        for (j = 0; j < max - min + 1; j++) {
            fsum += Data.allele_frequencies[i].freqs[j];
        }
        if (fsum > 1.01 || fsum < 0.99)
            fatalParseError(filename, msg_prefix,
                            "Allele frequencies do not sum to 1.", line);

        /* correct the small difference by altering the first allele frequency */
        Data.allele_frequencies[i].freqs[0] += 1.0 - fsum;
    }
}

void
DATAIOinit(void)
{
    int i;

    MAKE1DINT(Data.num_samples_in_population, Data.num_populations);
    MAKE1DINT(Data.num_ramets_in_population, Data.num_populations);
    MAKE1DINT(Data.locus_has_missing_data, Data.num_loci);
    MAKE1DDOUBLE(Data.TypingError, Data.num_loci);
    MAKE1DDOUBLE(Data._TypingErrorEst, Data.num_loci);
    MALLOC(Data.loci_ids, char*, Data.num_loci);

    if (Options.TypingErrorDefined >= 0)
        Data.TypingErrorAvg = Options.TypingErrorDefined;
    else
        Data.TypingErrorAvg = D_TYPINGERROR;

    for (i = 0; i < Data.num_loci; i++) {
        if (Options.TypingErrorDefined >= 0)
            Data.TypingError[i] = Options.TypingErrorDefined;
        else
            Data.TypingError[i] = D_TYPINGERROR;

        Data._TypingErrorEst[i] = Data.TypingError[i];
        MALLOC(Data.loci_ids[i], char, (D_PARSE_DESC_LENGTH+1));
        snprintf(Data.loci_ids[i], D_PARSE_DESC_LENGTH,"%i",(i+1));
    }
    MALLOC(Data.samples, SAMPLE *, Data.num_populations);
    MALLOC(Data.populations, POPULATION, Data.num_populations);

    if (Options.MinTyped < 0)
        Options.MinTyped = 1 + Data.num_loci / 2;

    if (Data.num_loci < Options.MinTyped)
        FATAL("Number of loci smaller than --mintyped.");
    Data.has_age_data = false;
    Data.has_missing_birth_data = false;
    Data.has_sex_data = false;
    Data.num_ramets = 0;
    Data.use_distances = false;
    Data.use_pedigree = false;
    Data._num_missing_data = 0;
    Data.sample_ids = NULL;
    Data.fs_matrix = NULL;
    Data.ignore_ary = NULL;
}

void
DATAIOcopyGenotypes(int ***source, int ***dest)
{
    int i, j, k;

    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < Data.num_loci; j++)
            for (k = 0; k < D_PLOIDY; k++)
                dest[i][j][k] = source[i][j][k];
}

int ***
DATAIOcloneGenotypes(void)
{
    int i, j, k, ***gcloned;

    MAKE3DINT(gcloned, Data.num_samples, Data.num_loci, 2, i, j);
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < Data.num_loci; j++)
            for (k = 0; k < D_PLOIDY; k++)
                gcloned[i][j][k] = Data.id_mapping[i].genotype_obs[j][k];
    return gcloned;
}

void
DATAIOdestroyGenotypes(int ***gcloned)
{
    int i, j;

    FREE3D(gcloned, Data.num_samples, Data.num_loci, i, j)
}

void
DATAIOcreateMapping(void)
{
    int i, j;

    /* now generate a mapping so that we can easily iterate over all 
     * samples in all populations */
    MALLOC(Data.id_mapping, SAMPLE, Data.num_samples);
    for (i = 0; i < Data.num_populations; i++) {
        for (j = 0; j < Data.num_samples_in_population[i]; j++) {
            Data.id_mapping[Data.samples[i][j].id] = Data.samples[i][j];
        }
    }
}

void
DATAIOdestroyMapping(void)
{
    FREE1D(Data.id_mapping);
}

void
DATAIOgenotypeout(char *filename)
{
    int i, j, k;
    FILE *out = NULL;
    char sex;

    FOPENW(out, filename);
    fprintf(out, "%d %d . %s\n", Data.num_populations, Data.num_loci,
            Data.description);
    for (i = 0; i < Data.num_populations; i++) {
        fprintf(out, "%d %s\n", Data.num_samples_in_population[i],
                Data.populations[i].description);
        for (j = 0; j < Data.num_samples_in_population[i]; j++) {
            fprintf(out, "%*s %d ", D_PARSE_DESC_LENGTH,
                    Data.samples[i][j].description,
                    Data.samples[i][j].ramets);

            if (Data.samples[i][j].birth < 0)
                fprintf(out, "? ");
            else
                fprintf(out, "%i ", Data.samples[i][j].birth);

            if (Data.samples[i][j].death < 0)
                fprintf(out, "? ");
            else
                fprintf(out, "%i ", Data.samples[i][j].death);

            sex = '?';
            if (Data.samples[i][j].sex == 1)
                sex = 'F';
            else if (Data.samples[i][j].sex == 2)
                sex = 'M';
            fprintf(out, "%c ", sex);

            for (k = 0; k < Data.num_loci; k++) {
                if (Data.samples[i][j].genotype_obs[k][0] > 0) {
                    if (Data.samples[i][j].genotype_obs[k][1] > 0)
                        fprintf(out, "%i.%i",
                                Data.samples[i][j].genotype_obs[k][0],
                                Data.samples[i][j].genotype_obs[k][1]);
                    else
                        fprintf(out, "%i.?",
                                Data.samples[i][j].genotype_obs[k][0]);
                } else {
                    if (Data.samples[i][j].genotype_obs[k][1] > 0)
                        fprintf(out, "?.%i",
                                Data.samples[i][j].genotype_obs[k][1]);
                    else
                        fprintf(out, "?.?");

                }
                if (k < Data.num_loci - 1)
                    fprintf(out, " ");
            }
            fprintf(out, "\n");
        }
    }
    FCLOSE(out);
}

void
DATAIOinfile(char *filename)
{
    int i, j, k, num_samples = 0, num_all_samples = 0, ng, res, line =
        0, a1, a2, tmp, num_loci_in_file;
    unsigned int l, l_old;
    FILE *in = NULL;
    char *p, c1, c2;
    char remaining[LINESIZE], input[LINESIZE], sep[1], msg[LINESIZE],
        msg_prefix[] = "INFILE; ";
    bool zero_allele = false, skip_missing;

    /* open file */
    if ((in = fopen(filename, "r")) == 0)
        FATAL("Cannot open input file!");

    /* read first line of the migrate formatted input file and parse it */
    FGETS;

    res =
        sscanf(input, "%d %d %1s %10s", &Data.num_populations, &Data.num_loci,
               sep, Data.description);

    if (Data.num_populations < 1 || Data.num_loci < 1)
        fatalParseError(filename, msg_prefix,
                        "Number sampling locations or loci too small.", line);
    if (Data.num_populations > D_PARSE_MAX_SAMPLING_LOCATIONS
        || Data.num_loci > D_PARSE_MAX_LOCI)
        fatalParseError(filename, msg_prefix,
                        "Number sampling locations or loci too high. Change global.h...",
                        line);

    /* we might change the true number of loci in the next lines, but we     *
     * need the true number for parsing loci ids                             */
    num_loci_in_file = Data.num_loci;

    if (Options.NumLoci > 0) {
        if (Options.NumLoci > Data.num_loci) {
            WARN("--numloci higher than the available number of loci. Ignoring this option.")
                Options.NumLoci = -1;
        } else {
            Data.num_loci = Options.NumLoci;
        }
    }

    if (res < 3)
        fatalParseError(filename, msg_prefix, "header not correct.", line);

    /* now we know with how many populations we are dealing with so          *
     * we can allocate space                                                 */
    DATAIOinit();

    /* now parse the data blocks (one for each population */
    for (i = 0; i < Data.num_populations; ++i) {
        FGETS;
        res = sscanf(input, "%d", &num_samples);
        if (!res) {
            if (i > 0)
            fatalParseError(filename, msg_prefix,
                            "number of individuals in sampling location not specified",
                            line);
            else {
                /* now assume the user provided loci ids */
                for (j=0; j< num_loci_in_file; j++) {
                    /* first already parsed */
                    if (j>0) {
                        FGETS;
                    }
                    /* --num_loci specified? then just skip the lines in     *
                     * the input file                                        */
                    if (j>= Data.num_loci) continue;

                    UTILScopyNonBlanks(input, Data.loci_ids[j], 0,
                                    D_PARSE_DESC_LENGTH);
                    UTILSstripNewline(Data.loci_ids[j], D_PARSE_DESC_LENGTH); 
                    UTILSremoveSpecialChars(Data.loci_ids[j],Data.loci_ids[j],D_PARSE_DESC_LENGTH);
                }
                FGETS;
                res = sscanf(input, "%d", &num_samples);
                if (res != 1) 
                    fatalParseError(filename, msg_prefix,
                                "number of individuals in sampling location not specified",
                                line);
            }        
        }    
        if (num_samples < 1)
            fatalParseError(filename, msg_prefix,
                            "number of individuals in sampling location too small.",
                            line);
        if (num_samples > D_PARSE_MAX_INDIVIDUALS)
            fatalParseError(filename, msg_prefix,
                            "number of individuals in sampling location too high. Change global.h if necessary...",
                            line);

        Data.populations[i].id = i;
        if ((p = strstr(input, " ")) != NULL) {
            strncpy(Data.populations[i].description, ++p,
                    D_PARSE_DESC_LENGTH);
            UTILSstripNewline(Data.populations[i].description,
                              D_PARSE_DESC_LENGTH);
        }
        Data.num_samples_in_population[i] = num_samples;
        Data.num_ramets_in_population[i] = 0;
        MALLOC(Data.samples[i], SAMPLE, num_samples);

        assert(Data.samples[i] != NULL);
        for (j = 0; j < num_samples; ++j) {
            MAKE2DINT(Data.samples[i][j].genotype_obs, Data.num_loci, 2, k);

            if (fgets(input, (int)sizeof(input), in) == NULL)
                fatalParseError(filename, msg_prefix,
                                "number of individuals/genotypes wrong.", line);
            line++;
            /* don't copy leading whitespaces */
            UTILScopyNonBlanks(input, Data.samples[i][j].description, 0,
                               D_PARSE_DESC_LENGTH);
            strcpy(remaining, input + D_PARSE_DESC_LENGTH);

            if ((p = strtok(remaining, " ")) == NULL) {
                if (line == 3)
                fatalParseError(filename, msg_prefix,
                 "Truncated line.\n\tMaybe invalid (numerical) loci ids?", line);
                else
                fatalParseError(filename, msg_prefix,
                                "Truncated line.", line);
            }    
            /* parse ramets */
            res = sscanf(p, "%d", &ng);
            if (res != 1) {
                if (line == 3)
                    fatalParseError(filename, msg_prefix,
                     "column <ramets> wrong.\n\tMaybe invalid (numerical) loci ids?", line);
                else 
                    fatalParseError(filename, msg_prefix, "column <ramets> wrong.", line);
            }    
            if (ng < 1 || ng > D_PARSE_MAX_RAMETS)
                fatalParseError(filename, msg_prefix,
                                "column <ramets> out of range. Change global.h...",
                                line);
            if ((p = strtok(NULL, " ")) == NULL)
                fatalParseError(filename, msg_prefix, "Could not parse column <ramets>.",line);
            Data.samples[i][j].ramets = ng;
            Data.num_ramets += ng;
            Data.num_ramets_in_population[i] += ng;
            /* parse age */
            res = sscanf(p, "%d", &Data.samples[i][j].birth);
            if (res != 1) {
                res = sscanf(p, "%c", &c1);
                if (res != 1 || c1 != '?')
                    fatalParseError(filename, msg_prefix,
                                    "Character in column <birth> wrong. Must be 0-9 or ?.",
                                    line);
                /* c1 is ?, unknown, so set age to negative value */
                Data.samples[i][j].birth = -1;
                Data.has_missing_birth_data = true;
            } else if (Options.IgnoreAge) {
                Data.samples[i][j].birth = -1;
                Data.has_missing_birth_data = true;
            } else {
                Data.has_age_data = true;
            }
            if ((p = strtok(NULL, " ")) ==NULL)
                    fatalParseError(filename, msg_prefix,
                                    "Could not parse column <death>.",line);
            res = sscanf(p, "%d", &Data.samples[i][j].death);
            if (res != 1) {
                res = sscanf(p, "%c", &c1);
                if (res != 1 || c1 != '?')
                    fatalParseError(filename, msg_prefix,
                                    "Character in column <death> wrong. Must be 0-9 or ?.",
                                    line);
                /* c1 is ?, unknown, so set age to negative value */
                Data.samples[i][j].death = -1;
            } else if (Options.IgnoreAge) {
                Data.samples[i][j].death = -1;
            } else {
                Data.has_age_data = true;
            }
            p = strtok(NULL, " ");
            /* parse sex */
            res = sscanf(p, "%c", &c2);
            if (res != 1)
                fatalParseError(filename, msg_prefix,
                                "Character in column <sex> wrong. Must be M,F or ?", line);
            if (Options.IgnoreSex) {
                Data.samples[i][j].sex = -1;
            }    
            else {    
                switch (toupper(c2)) {
                    case 'F':
                        Data.samples[i][j].sex = 1;
                        Data.has_sex_data = true;
                        break;
                    case 'M':
                        Data.samples[i][j].sex = 2;
                        Data.has_sex_data = true;
                        break;
                    case '?':
                        Data.samples[i][j].sex = -1;
                        break;
                    default:
                        fatalParseError(filename, msg_prefix,
                                        "sex must be M, F or ?", line);
                }
            }
            p = strtok(NULL, " ");

            Data.samples[i][j].id = num_all_samples++;
            if (!addSampleIDtoHash
                (Data.samples[i][j].id, Data.samples[i][j].description)) {
                snprintf(msg, sizeof(msg), "id (%s) not unique",
                         Data.samples[i][j].description);
                fatalParseError(filename, msg_prefix, msg, line);
            }

            Data.samples[i][j].population = &Data.populations[i];
            Data.samples[i][j].typed_loci = 0;

            /* parse marker data */
            for (k = 0; k < Data.num_loci; ++k) {
                if (p == NUL) {
                    snprintf(msg, sizeof(msg), "Expecting %i loci, but got only %i", Data.num_loci, k);
                    fatalParseError(filename, msg_prefix, msg, line);
                }    
                res = sscanf(p, "%d%c%d", &a1, sep, &a2);
                if (res == 2) { /* ok, parsing of the second allele failed */
                    /* second allele a question mark (=missing data)?       */
                    res = sscanf(p, "%d%c%c", &a1, sep, &c2);
                    if (res == 3 && c2 == '?')
                        a2 = -1;
                    else {
                        snprintf(msg, sizeof(msg), "Second allele of locus %i not a number or a '?'", k);
                        fatalParseError(filename, msg_prefix, msg,line);
                    }    
                } /* parsing of both alleles failed */
                else if (res == 0) {
                    /* first allele missing? */
                    res = sscanf(p, "%c%c%d", &c1, sep, &a2);
                    if (res == 3) {
                        a1 = -1;
                        if (c1 != '?') {
                            snprintf(msg, sizeof(msg), "First allele of locus %i not a number or a '?'", k);
                            fatalParseError(filename, msg_prefix, msg,line);
                        }
                    } else {
                        /* both missing? */
                        res = sscanf(p, "%c%c%c", &c1, sep, &c2);
                        if (res == 3) {
                            a1 = -1;
                            a2 = -1;
                            if (c1 != '?' || c2 != '?') {
                                snprintf(msg, sizeof(msg), "Alleles of locus %i not numbers or '?'s", k);
                                fatalParseError(filename, msg_prefix,msg, line);
                            }    
                        } else {
                            snprintf(msg, sizeof(msg), "Given up while parsing locus %i", k);
                            fatalParseError(filename, msg_prefix,msg,line);
                        }    
                    }
                }
                if (a2 < a1)
                    SWAP(a1, a2, tmp);

                if (a1 < 0 || a2 < 0) {
                    Data.locus_has_missing_data[k] = 1;
                    Data._num_missing_data++;
                    assert(a2 >= a1);
                    if (a2 < 0)
                        Data._num_missing_data++;
                } else
                    Data.samples[i][j].typed_loci++;

                if (a1 == 0 || a2 == 0) zero_allele = true;

                if (a1 > D_PARSE_MAX_ALLELE_VALUE || a2 > D_PARSE_MAX_ALLELE_VALUE) {
                    snprintf(msg, sizeof(msg), "Value of allele at locus %i too big!", k);
                    fatalParseError(filename, msg_prefix,msg, line);
                }                
                Data.samples[i][j].genotype_obs[k][0] = a1;
                Data.samples[i][j].genotype_obs[k][1] = a2;
                p = strtok(NULL, " ");
            }
        }
    }
    /* close file handle */
    FCLOSE(in);

    Data.num_samples = num_all_samples;

    MAKE2DBOOL(Data.ignore_ary, Data.num_samples, Data.num_loci, i);
    DATAIOcreateMapping();
    
    if (!Data.has_sex_data && !Options.MonoeciousDefined) Options.Monoecious = true;
    if (Data.has_sex_data && Options.Monoecious) 
        WARN("Sex data specified and monoecious options set.");

    if (Data._num_missing_data > 0) {
        l = 0;
        MALLOC(Data._missing_data, LOCUS_COORD, Data._num_missing_data);
        for (i = 0; i < Data.num_samples; i++) {
            if (Data.id_mapping[i].typed_loci < Options.MinTyped) 
                skip_missing = true;
            else
                skip_missing = false;
            for (j = 0; j < Data.num_loci; j++) {
                l_old = l;
                for (k = 0; k < D_PLOIDY; k++)
                    if (Data.id_mapping[i].genotype_obs[j][k] < 0) {
                        assert(l < Data._num_missing_data);
                        if (skip_missing) {
                            Data._num_missing_data--;
                        }
                        else {
                            Data._missing_data[l].id = i;
                            Data._missing_data[l].locus = j;
                            Data._missing_data[l].allele = k;
                            Data._missing_data[l].both = false;
                            l++;
                        }
                    }
                if (l - l_old == 2) {
                    Data._missing_data[l-2].both = true;
                    Data._missing_data[l-1].both = true;
                }        
            }
        } 
    }

    /* calculate log(n!) for 0..2n */
    UTILSfactlnInit(Data.num_samples * 2);

    if (Data.has_age_data && Options.FemRepro.min == D_FEMREPROMIN
        && Options.FemRepro.max == D_FEMREPROMAX)
        WARN("Data has age data and --femrepro is not set.");

    if (Data.has_age_data && Options.MaleRepro.min == D_MALEREPROMIN &&
        Options.MaleRepro.max == D_MALEREPROMAX)
        WARN("Data has age data and --malerepro is not set.");

    if (zero_allele)
        WARN("There are alleles coded as '0' in the data. Do you mean '?' (missing value) instead?");

}

static void
calcMismatches(int **child, int **parent, int *mismatching, int *mismatch_at_locus)
{
    int i;
    double p = 0.;

    *mismatching = 0;

    for (i = 0; i < Data.num_loci; i++) {

        if (child[i][0] < 0 && child[i][1] < 0) continue;

        /* all typed? then count as known relationship */
        if (child[i][0] >= 0 && child[i][1] >= 0 && parent[i][0] >= 0
            && parent[i][1] >= 0)
            Data.in_pedigree_relationships[i]++;

        p = LODcalcTransProbDyad(child[i], parent[i], i);
        if (CMP_DBL(p,0.)) {
            (*mismatching)++;
            mismatch_at_locus[i] = 1;
        }
    }
}

void
DATAIOestimateTypingError(void)
{
    int i, j, k, mismatches, max_mm = 0, relationships = 0;
    double typingerror = D_TYPINGERROR, sum_te = 0., min_relationships;

    /* if we don't have a pedigree we can't estimate the typing error and if
     * the user aleady defined an error rate, we don't want to overwrite this
     * value */
    if (!Options.MaxMismatchingDefined) {
        if (Options.TypingErrorDefined >= 0)
            typingerror = Options.TypingErrorDefined;
        Options.MaxMismatchingDyad = (short)ceil(Data.num_loci * 4 * typingerror);
        Options.MaxMismatchingTriple = Options.MaxMismatchingDyad + 1;
    }

    /* if there are no relationships known, we cannot look for 
     * mismatches */
    if (!Data.use_pedigree)
        return;

    for (i = 0; i < Data.num_loci; i++)
        Data.in_pedigree_relationships[i] = 0;

    for (i = 0; i < Data.num_samples; i++) {
        for (j = 0; j < Data.num_loci; j++) {
            Data.in_pedigree_mismatches[i][0][j] = 0;
            Data.in_pedigree_mismatches[i][1][j] = 0;
        }
        for (j = 0; j < 2; j++) {
            k = Data.in_pedigree[i][j];
            if (k >= 0) {
                calcMismatches(Data.id_mapping[i].genotype_obs,
                               Data.id_mapping[k].genotype_obs, &mismatches,
                               Data.in_pedigree_mismatches[i][j]);
                if (mismatches > max_mm)
                    max_mm = mismatches;
            }
        }
    }

    if (Options.TypingErrorDefined >= 0)
        return;

    for (i = 0; i < Data.num_loci; i++) {
        relationships += Data.in_pedigree_relationships[i];
        mismatches = 0;
        for (j = 0; j < Data.num_samples; j++) {
            mismatches += Data.in_pedigree_mismatches[j][0][i];
            mismatches += Data.in_pedigree_mismatches[j][1][i];
        }
        /* see marshall 1998, appendix */
        Data.TypingError[i] =
            (1. / (double)(2. * Data.exclp_sp[0][i])) * (mismatches /
                                                         (double)Data.in_pedigree_relationships[i]);
        Data._TypingErrorEst[i] = Data.TypingError[i];
        sum_te += Data.TypingError[i];
        /* always allow a small error, use should turn this off when necessary */
        if (isnan(Data.TypingError[i])
            || Data.TypingError[i] <= D_TYPINGERRORMIN)
            Data.TypingError[i] = D_TYPINGERRORMIN;
    }
    Data.TypingErrorAvg = sum_te / (double)Data.num_loci;
    min_relationships = (D_ERRORZ / D_ERRORE) * (D_ERRORZ / D_ERRORE);

    for (i = 0; i < Data.num_loci; i++)
        /* we don't have enough relationships? */
        if (Data.in_pedigree_relationships[i] < (int)ceil(min_relationships)) {
            /* then assume that error rate is constant across loci */
            if (relationships > (int)ceil(min_relationships / Data.num_loci)) {
                Data.TypingError[i] =
                    MAX(Data.TypingErrorAvg, D_TYPINGERRORMIN);
            } else {
                /* else use default error rate                         */
                Data.TypingError[i] = D_TYPINGERROR;
            }
        }

    if (!Options.MaxMismatchingDefined) {
        Options.MaxMismatchingTriple =
            MAX(max_mm, Options.MaxMismatchingTriple);
        Options.MaxMismatchingDyad = MAX(max_mm, Options.MaxMismatchingDyad);
    }
}

/* returns 1 if the genotypes are identical */
static inline bool
compareGenotypes(int id1, int id2)
{
    int i, j;

    assert(id1 != id2);
    for (i = 0; i < Data.num_loci; i++) {
        /* assumes that alleles are sorted, DATIOinfile does this */
        for (j = 0; j < D_PLOIDY; j++)
            if (Data.id_mapping[id1].genotype_obs[i][j] !=
                Data.id_mapping[id2].genotype_obs[i][j])
                return false;
    }
    return true;
}

void
DATAIOcheckData(FILE * fp)
{
    int i, j;

    bool *dup, dups_found = false, dups_found_ind;
    MAKE1DBOOL(dup, Data.num_samples);

    for (i = 0; i < Data.num_samples; i++) 
        dup[i] = false;

    for (i = 0; i < Data.num_samples; i++) {
        dups_found_ind = false;
        /* already identified as duplicate?*/
        if (dup[i]) continue;

        for (j = i+1; j < Data.num_samples; j++) {
            if (compareGenotypes(i, j)) {
                dup[j] = true;
                if (!dups_found) 
                    fprintf(fp, "*** Identical genotypes ***\n\n");
                
                if (!dups_found_ind)
                    fprintf(fp, "%-*s", D_PARSE_DESC_LENGTH,
                            Data.id_mapping[i].description);
                dups_found_ind = true;
                dups_found = true;

                fprintf(fp, "  %-*s", D_PARSE_DESC_LENGTH,
                        Data.id_mapping[j].description);
            }
        }
        if (dups_found_ind)
            fprintf(fp, "\n");
    }
    FREE1D(dup);
    if (dups_found)
        fprintf(fp, "\n\n");

}

void DATAIOdumpMismatches(void)
{
    int i, j, k, l, l1_het, l2_het, mm_complete, mm_het_het,
        mm_hom_hom, mm_het_hom, mm_found = 0;
    bool found;
    unsigned int sum, total = 0;
    double percentage;
    double sum_d;
    FILE *out;
    FOPENW(out, Options.MismatchOutfilename);
    DATAIOdumpDatasetDetails(out);
    if (Data.use_pedigree) {
        for (i = 0; i < Data.num_samples; i++) {
            for (j = 0; j < 2; j++)
                if (Data.in_pedigree[i][j] >= 0) {
                    k = Data.in_pedigree[i][j];
                    for (l = 0; l < Data.num_loci; l++)
                        if (Data.in_pedigree_mismatches[i][j][l] > 0) {
                            if (mm_found == 0) {

                                fprintf(out, "*** Known parent-offspring mismatches ***\n\n");
                                fprintf(out, "%-*s  %-12s  %-13s  %-15s  %-13s\n",D_PARSE_DESC_LENGTH, "Locus", "Offspring ID", "Genotype", "Known parent ID", "Genotype");
                            }
                            fprintf(out, "%-*s  %-12s  %-6i %-6i  %-15s  %-6i %-6i\n", D_PARSE_DESC_LENGTH, Data.loci_ids[l], Data.id_mapping[i].description,
                                    Data.id_mapping[i].genotype_obs[l][0],
                                    Data.id_mapping[i].genotype_obs[l][1],Data.id_mapping[k].description,
                                    Data.id_mapping[k].genotype_obs[l][0],
                                    Data.id_mapping[k].genotype_obs[l][1]);
                            mm_found++;
                        }
                }
        }
        if (mm_found > 0)
            fprintf(out, "\nTOTAL: %i\n\n\n",mm_found);
    
        fprintf(out, "*** Summary Observed Mismatches ***\n\n");
        fprintf(out,
                "%-*s  Hom/Hom  Het/Het  Hom/Het   Sum   Relationships  Detection Prob.  Est. Error\n", D_PARSE_DESC_LENGTH, "Locus");
        for (i = 0; i < Data.num_loci; i++) {
            mm_complete = 0;
            mm_het_hom = 0;
            mm_het_het = 0;
            mm_hom_hom = 0;
            for (k = 0; k < 2; k++)
                for (j = 0; j < Data.num_samples; j++) {
                    l = Data.in_pedigree[j][k];
                    if (l < 0 || !Data.in_pedigree_mismatches[j][k][i])
                        continue;
                    mm_complete++;
                    l1_het = 1;
                    l2_het = 1;
                    if (Data.id_mapping[j].genotype_obs[i][0] ==
                        Data.id_mapping[j].genotype_obs[i][1])
                        l1_het = 0;
                    if (Data.id_mapping[l].genotype_obs[i][0] ==
                        Data.id_mapping[l].genotype_obs[i][1])
                        l2_het = 0;
                    if (l1_het && l2_het)
                        mm_het_het++;
                    else if (!l1_het && !l2_het)
                        mm_hom_hom++;
                    else
                        mm_het_hom++;
                    assert(mm_complete ==
                           (mm_het_het + mm_het_hom + mm_hom_hom));
                }
            fprintf(out, "%-*s  %7i  %7i  %7i  %5i  %13i %15.4f  %10.4f\n", D_PARSE_DESC_LENGTH, Data.loci_ids[i],
                    mm_hom_hom, mm_het_het, mm_het_hom, mm_complete,
                    Data.in_pedigree_relationships[i], Data.exclp_sp[0][i], Data._TypingErrorEst[i]);
        }
        fprintf(out, "\n\n");
    }
    fprintf(out, "*** MCMC Mismatches ***\n\n");
    fprintf(out, "Only mismatches with a percentage >= %.2f are listed.\n\n", D_MISMATCHOUTTH);

    for (i = 0; i < Data.num_loci; i++) {
        found = false;
        total = 0;
        for (j=0; j<Data.num_samples; j++) {
        
            if (Probs.num_candidates[j] == 0) continue;

            sum = 0; sum_d = 0;
            if (Probs.mh_sampled_pedigrees > 0) {
                for (k=0;k< Probs.num_posteriors[j]; k++) 
                    if (LODhasMismatchLocus(i,Data.id_mapping[j].genotype_obs,Probs.posteriors[j][k].v, Probs.posteriors[j][k].w, Data.ignore_ary[j])) 
                        sum += Probs.posteriors[j][k].observed;
                percentage = sum/(double)(Probs.mh_sampled_pedigrees);
            }
            else {
                for (k=0;k< Probs.num_posteriors[j]; k++) 
                    if (LODhasMismatchLocus(i,Data.id_mapping[j].genotype_obs,Probs.posteriors[j][k].v, Probs.posteriors[j][k].w, Data.ignore_ary[j])) 
                       sum_d += exp(Probs.posteriors[j][k].p_opt);
                percentage = sum_d;
            }    

            if (percentage < D_MISMATCHOUTTH) continue;

            total++;

            if (!found) {
                fprintf(out, "*** Locus %s ***\n\n", Data.loci_ids[i]);
                found = true;
            }
            fprintf(out, "%-10s %10.4f %%\n", Data.id_mapping[j].description,percentage);
        }
        if (found) {
            fprintf(out, "\nTOTAL: %u\n\n", total);
            FFLUSH(out);
        }
    }
    FCLOSE(out);
}

void
DATAIOdump(FILE * fp)
{
    int i, j;
    SAMPLE s;

    fprintf(fp, "%i samples in %i populations.\n", Data.num_samples,
            Data.num_populations);
    for (i = 0; i < Data.num_samples; i++) {
        s = Data.id_mapping[i];
        fprintf(fp, "%4d %i ", i, s.ramets);
        for (j = 0; j < Data.num_loci; j++) {
            fprintf(fp, "%i.%i ", s.genotype_obs[j][0], s.genotype_obs[j][1]);
        }
        fprintf(fp, "\n");
    }
}

void
DATAIOdumpDatasetDetails(FILE * fp)
{
    time_t now;

    (void)time(&now);
    fprintf(fp, "Software Version            : %s %s\n", PACKAGE,
            PACKAGE_VERSION);
    fprintf(fp, "Dataset                     : %s\n", Data.description);
    fprintf(fp, "Time of Data Analysis       : %s\n", ctime(&now));
    fprintf(fp, "\n");
}

GraphvizCluster
DATAIOgraphvizCluster(void)
{
    int i, j;
    GraphvizCluster cluster;

    cluster.num_cluster = Data.num_populations;
    MALLOC(cluster.num_vertices_in_cluster, int, Data.num_populations);
    MALLOC(cluster.ids, int *, Data.num_populations);
    MALLOC(cluster.labels, char *, Data.num_populations);

    for (i = 0; i < Data.num_populations; i++) {
        MALLOC(cluster.ids[i], int, Data.num_samples_in_population[i]);

        cluster.num_vertices_in_cluster[i] =
            Data.num_samples_in_population[i];
        for (j = 0; j < Data.num_samples_in_population[i]; j++) {
            cluster.ids[i][j] = Data.samples[i][j].id;
        }
        MALLOC(cluster.labels[i], char, 51);

        strncpy(cluster.labels[i], Data.samples[i][0].population->description,
                50);
    }
    return cluster;
}

void
DATAIOgraphvizClusterDestroy(GraphvizCluster cluster)
{
    int i;

    FREE2D(cluster.ids, Data.num_populations, i);
    FREE2D(cluster.labels, Data.num_populations, i);
    FREE(cluster.num_vertices_in_cluster);
}

static void
hlSL(FILE *fp)
{
    int i,j;
    fprintf(fp, "\n ");
    for (i = 0; i < ( Data.num_populations + 1); i++) {
        fprintf(fp, "+");
        for (j = 0; j < (D_PARSE_DESC_LENGTH + 2); j++) 
            fprintf(fp, "-");
    }

    fprintf(fp, "+");
}

void
DATAIOdumpSamplingLocations(FILE * fp)
{
    int i,j;

    if (!Data.use_distances)
        return;

    fprintf(fp, "*** Sampling locations ***\n");

    hlSL(fp);
    fprintf(fp, "\n | %-*s |", D_PARSE_DESC_LENGTH, " ");
    for (i = 0; i < Data.num_populations; i++) 
        fprintf(fp, " %-*s |", D_PARSE_DESC_LENGTH, Data.populations[i].description);

    hlSL(fp);
    for (i = 0; i < Data.num_populations; i++) {
        fprintf(fp,"\n | %-*s |", D_PARSE_DESC_LENGTH, Data.populations[i].description);
        for (j = 0; j < Data.num_populations; j++) 
            fprintf(fp, " %*.2f |", D_PARSE_DESC_LENGTH, Data.population_distances[i][j]);
    }
    hlSL(fp);
    fprintf(fp, "\n\n\n");
}

void
DATAIOaddGVNodeStyle(Dag D)
{
    int i;
    GraphvizVertex gv;

    for (i = 0; i < Data.num_samples; i++) {
        DAGlabelV(D, i, Data.id_mapping[i].description, D_PARSE_DESC_LENGTH);
        if (Data.id_mapping[i].sex == 1)
            strcpy(gv.shape, "circle");
        else if (Data.id_mapping[i].sex == 2)
            strcpy(gv.shape, "box");
        else
            strcpy(gv.shape, "diamond");
        DAGgraphvizV(D, i, gv);
    }
}

void
DATAIOdumpGenotypesCERVUS(void)
{
    int i, j, k;
    FILE *fp;

    FOPENW(fp, Options.GenotypeOutfilenameCERVUS);
    fprintf(fp, "Individual ID,Sex");

    for (i = 0; i < Data.num_loci; i++)
        fprintf(fp, ",%sa,%sb", Data.loci_ids[i], Data.loci_ids[i]);

    for (i = 0; i < Data.num_samples; i++) {
        fprintf(fp, "\n%s", Data.id_mapping[i].description);
        if (Data.id_mapping[i].sex < 0)
            fprintf(fp, ",");
        else
            fprintf(fp, ",%i", Data.id_mapping[i].sex);
        for (j = 0; j < Data.num_loci; j++) {
            for (k = 0; k < D_PLOIDY; k++)
                if (Data.id_mapping[i].genotype_obs[j][k] < 0)
                    fprintf(fp, ",");
                else
                    fprintf(fp, ",%i", Data.id_mapping[i].genotype_obs[j][k]);
        }
    }
    FCLOSE(fp);
}

void
DATAIOdumpGenotypesGenepop(void)
{
    int i, j, k, l,*allele_strlens;
    FILE *fp;
    char allele[D_PARSE_MAX_ALLELE_VALUE_STRLEN+1];

    FOPENW(fp, Options.GenotypeOutfilenameGenepop);
    if (strlen(Data.description) > 0)
        fprintf(fp, "%s\n",Data.description);
    else
        fprintf(fp, "Title\n");
    
    MAKE1DINT(allele_strlens, Data.num_loci);

    for (i = 0; i < Data.num_loci; i++) {
        fprintf(fp, "%s\n", Data.loci_ids[i]);
        /* get length of allele with highest value (2,3 or 4?) */
        snprintf(allele, D_PARSE_MAX_ALLELE_VALUE_STRLEN+1, "%d", Data.allele_frequencies[i].max);
        allele_strlens[i] = MAX(2,strlen(allele));
    }

    for (i = 0; i < Data.num_populations; i++) {
        fprintf(fp, "Pop\n");
        for (j = 0; j < Data.num_samples_in_population[i]; j++) {
            if (Data.num_populations > 1 && strlen(Data.populations[i].description) > 0)
                fprintf(fp, "%*s %*s , ", D_PARSE_DESC_LENGTH, Data.populations[i].description, D_PARSE_DESC_LENGTH, Data.samples[i][j].description);
            else
                fprintf(fp, "%*s , ", D_PARSE_DESC_LENGTH, Data.samples[i][j].description);
            
            for (k = 0; k < Data.num_loci; k++) { 
                fprintf(fp, " ");
                for (l = 0; l < D_PLOIDY; l++) 
                    fprintf(fp, "%0*i", allele_strlens[k], MAX(Data.samples[i][j].genotype_obs[k][l], 0));
            }
            fprintf(fp, "\n");
        }
    }
    FCLOSE(fp);
    FREE1D(allele_strlens);
}

void 
DATAIOdumpGenotypesRMES(void)
{
    int i, j, a;
    FILE *fp;

    FOPENW(fp, Options.GenotypeOutfilenameRMES);

    fprintf(fp, "%d\r\n",Data.num_populations);
    
    for (i = 0; i < Data.num_populations; i++) 
        fprintf(fp, "%*s\r\n%d\r\n%d\r\n", D_PARSE_DESC_LENGTH, Data.populations[i].description, Data.num_samples_in_population[i],Data.num_loci);

    for (i = 0; i < Data.num_samples; i++) {
        for (j=0; j<Data.num_loci;j++) {
            /* missing values coded as 99 as recommended in the manual,  *
             * het = 1, hom = 0                                          */
            if (Data.id_mapping[i].genotype_obs[j][0] < 0 || Data.id_mapping[i].genotype_obs[j][1] < 0) a = -99;
            else if (Data.id_mapping[i].genotype_obs[j][0] == Data.id_mapping[i].genotype_obs[j][1]) a = 0;
            else a=1;

            fprintf(fp, "%d ", a);
        }
        fprintf(fp, "\r\n");
    }
    FCLOSE(fp);
}

void 
DATAIOdumpGenotypesParente(void)
{
    int i, j, k;
    FILE *fp;
    char sex;

    FOPENW(fp, Options.GenotypeOutfilenamePARENTE);
    
    fprintf(fp, "\tSex\tYear of birth\tYear of death\tMother");

    for (i=0;i<Data.num_loci;i++)
        fprintf(fp, "\t%s1\t%s2",Data.loci_ids[i],Data.loci_ids[i]);

    fprintf(fp, "\r\n");

    for (i = 0; i < Data.num_samples; i++) {
        sex = '?';
        if (Data.id_mapping[i].sex == 1)
            sex = 'F';
        else if (Data.id_mapping[i].sex == 2)
            sex = 'M';

        fprintf(fp, "%*s\t%c", D_PARSE_DESC_LENGTH, Data.id_mapping[i].description, sex);
        if (Data.id_mapping[i].birth < 0)
            fprintf(fp, "\t%c", '?');
        else
            fprintf(fp, "\t%i", Data.id_mapping[i].birth);

        if (Data.id_mapping[i].death < 0)
            fprintf(fp, "\t%c", '?');
        else
            fprintf(fp, "\t%i", Data.id_mapping[i].death);
        
        if (HAS_KNOWN_MOTHER(i))
            fprintf(fp, "\t%*s", D_PARSE_DESC_LENGTH,Data.id_mapping[Data.in_pedigree[i][0]].description);
        else 
            fprintf(fp, "\t?");

        for (j=0; j<Data.num_loci;j++) {
            for (k=0;k<D_PLOIDY;k++)
                if (Data.id_mapping[i].genotype_obs[j][k] < 0) 
                    fprintf(fp,"\t?");
                else 
                    fprintf(fp, "\t%i",Data.id_mapping[i].genotype_obs[j][k]);
        }
        fprintf(fp, "\r\n");
    }
    FCLOSE(fp);
}

void
DATAIOdumpOffspringCERVUS(void)
{
    int i, j;
    FILE *fp;

    FOPENW(fp, Options.OffspringOutfilenameCERVUS);
    fprintf(fp,
            "Offspring ID,Year of birth,Known mother ID,Candidate father IDs");

    for (i = 0; i < Data.num_samples; i++) {
        if (Probs.num_candidates[i] <= 0)
            continue;

        fprintf(fp, "\n%s", Data.id_mapping[i].description);
        if (Data.id_mapping[i].birth < 0)
            fprintf(fp, ",");
        else
            fprintf(fp, ",%i", Data.id_mapping[i].birth);

        if (HAS_KNOWN_MOTHER(i))
            fprintf(fp, ",%s",
                    Data.id_mapping[Data.in_pedigree[i][0]].description);
        else
            fprintf(fp, ",");
        
        if (Probs.num_candidates[i] == Data.num_samples - 1)
            for (j = 0; j < Data.num_samples; j++) {
                if (i == j) continue;
                fprintf(fp, ",%s",
                        Data.id_mapping[j].description);
            }
        else
            for (j = 0; j < Probs.num_candidates[i]; j++)
                fprintf(fp, ",%s",
                        Data.id_mapping[Probs.candidates[i][j]].description);
    }
    FCLOSE(fp);
}

void
DATAIOinitFullsibMatrix(void)
{
    int i, j, t;

    /* already initialized? */
    if (Data.fs_matrix != NULL)
        return;

    MAKETRIANGULAR(Data.fs_matrix, FULLSIB, Data.num_samples);

    for (i = 0; i < Data.num_samples; i++) {
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            Data.fs_matrix[t].fullsibs = false;
            Data.fs_matrix[t].halfsibs = false;
            Data.fs_matrix[t].indirect = false;
        }
    }
}

void
DATAIOdestroyPreMCMC(void)
{
    if (Data.fs_matrix != NULL) {
        FREETRIANGULAR(Data.fs_matrix);
    }
    UTILSfactlnDestroy();
}

void
DATAIOdestroyPostMCMC(void)
{
    int i, j, k;
    DESCRIPTION_HASH *s;

    for (i = 0; i != Data.num_populations; ++i) {
        for (j = 0; j != Data.num_samples_in_population[i]; ++j) {
            if (Data.sample_ids != NULL) {
                HASH_FIND_STR(Data.sample_ids, Data.samples[i][j].description,
                              s);
                HASH_DEL(Data.sample_ids, s);
                FREE(s);
            }
            FREE2D(Data.samples[i][j].genotype_obs, Data.num_loci, k);
            Data.samples[i][j].population = NULL;
        }
    }
    FREE1D(Data._TypingErrorEst);
    FREE1D(Data.populations);
    FREE1D(Data.locus_has_missing_data);
    FREE2D(Data.samples, Data.num_populations, i);
    if (Data.use_distances)
        FREE2D(Data.population_distances, Data.num_populations, i);
    if (Data.use_pedigree) {
        FREE1D(Data.in_pedigree_relationships);
        FREE2D(Data.in_pedigree, Data.num_samples, i);
        FREE3D(Data.in_pedigree_mismatches, Data.num_samples, 2, i, j);
    }
    DATAIOdestroyMapping();
    FREE1D(Data.num_samples_in_population);
    FREE1D(Data.num_ramets_in_population);
    FREE1D(Data.TypingError);
    FREE2D(Data.loci_ids,Data.num_loci, i);
    if (Data.ignore_ary != NULL)
        FREE2D(Data.ignore_ary,Data.num_samples, i);
    if (Data._num_missing_data > 0) {
        FREE1D(Data._missing_data);
    }    
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
