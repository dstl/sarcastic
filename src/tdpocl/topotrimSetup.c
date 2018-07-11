/***************************************************************************
 *
 *       Module:    topotrimSetup.c
 *      Program:    tdp2
 *   Created by:    Darren Muff on 14/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/

#include <stdio.h>
#include "tdp.h"

void
topotrim_setup(FILE **topotrim_fp, char **input_filename, char **surface_fname,
               int *NX, int *NY, int *xstart, int *ystart)
/* function to parse the Topotrim file and get the key parameters from its header */
{
    
    char * line_ref, *line_mis, *line_surf;
    char * topotrim_filename;
    char * ref_fname, * mis_fname;
    int process_reference ;
    
    topotrim_filename = input_string("Input TTM data", "InputTTMName", "The name of the output file from TopoTrim",
                                     "first_second.ttm");
    
    /* Lets open the file and check that the first 3 lines start the way they should */
    *topotrim_fp = fopen(topotrim_filename, "r");
    if (!*topotrim_fp) {
        fprintf(stderr, "Failed to open file %s for read!\n", topotrim_filename);
        perror("File open fail");
        exit(12345);
    }
    
    char *refmisstr;
    do {
        refmisstr = input_string("Which image to process? (r)eference or (m)ission ?", "ProcessRef",
                                 "enter \'r\' to process the reference image or \'m\' for the misison image","r");
    }while (refmisstr[0] != 'r' && refmisstr[0] != 'R' && refmisstr[0] != 'm' && refmisstr[0] != 'M');
    
    if(refmisstr[0] == 'r' || refmisstr[0] == 'R'){
        process_reference = 1 ;
    }else{
        process_reference = 0;
    }
    
    line_ref = calloc(TTM_LINE_SIZE, 1);
    fgets(line_ref, TTM_LINE_SIZE, *topotrim_fp);
    if (strncmp(line_ref, "Reference CPHD File is :", 24) != 0) {
        fprintf(stderr, "Badly formatted TTM file!\n Expected line to start with \"Reference CPHD File is :\", actually was \"%s\"\n\n", line_ref);
        exit(12346);
    }
    ref_fname = strstr(line_ref, ":");
    if (!ref_fname) {
        fprintf(stderr, "Weird - couldn't find the \":\" in line_ref (%s)\n", line_ref);
        exit(12347);
    }
    ref_fname++;   //  ref_name[0] now == ':', so move it along by 1 to point to the actual string we want
    ref_fname[strlen(ref_fname)-1] = 0; // Splat the \n at the end of the string
    
    line_mis = calloc(TTM_LINE_SIZE, 1);
    fgets(line_mis, TTM_LINE_SIZE, *topotrim_fp);
    if (strncmp(line_mis, "Mission CPHD File is   :", 24) != 0) {
        fprintf(stderr, "Badly formatted TTM file!\n Expected line to start with \"Mission CPHD File is   :\", actually was \"%s\"\n\n", line_mis);
        exit(12348);
    }
    mis_fname = strstr(line_mis, ":");
    if (!mis_fname) {
        fprintf(stderr, "Weird - couldn't find the \":\" in line_mis (%s)\n", line_mis);
        exit(12349);
    }
    mis_fname++;
    mis_fname[strlen(mis_fname)-1] = 0; // Splat the \n at the end of the string
    
    line_surf = calloc(TTM_LINE_SIZE, 1);
    fgets(line_surf, TTM_LINE_SIZE, *topotrim_fp);
    if (strncmp(line_surf, "Surface file is        :", 24) != 0) {
        fprintf(stderr, "Badly formatted TTM file!\n Expected line to start with \"Surface file is        :\", actually was \"%s\"\n\n", line_surf);
        exit(12350);
    }
    *surface_fname = strdup(strstr(line_surf, ":"));
    if (!(*surface_fname)) {
        fprintf(stderr, "Weird - couldn't find the \":\" in line_surf (%s)\n", line_surf);
        exit(12351);
    }
    (*surface_fname)++;
    (*surface_fname)[strlen(*surface_fname)-1] = 0; // Splat the \n at the end of the string
    
    printf("\nThe reference filename: %s\n", ref_fname);
    printf("The mission filename: %s\n", mis_fname);
    printf("The surface filename: %s\n\n", *surface_fname);
    
    if (process_reference) {
        *input_filename = strdup(ref_fname);
    } else {
        *input_filename = strdup(mis_fname);
    }
    
    int t1,t2,t3,t4;
    fread(&t1, sizeof(int), 1, *topotrim_fp) ;
    fread(&t2, sizeof(int), 1, *topotrim_fp) ;
    fread(&t3, sizeof(int), 1, *topotrim_fp) ;
    fread(&t4, sizeof(int), 1, *topotrim_fp) ;
    *NX = t1;
    *NY = t2;
    *xstart = t3;
    *ystart = t4;
    
    free(line_mis);
    free(line_ref);
    free(line_surf);
    
    
}
