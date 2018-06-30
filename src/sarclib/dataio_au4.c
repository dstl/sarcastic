/***************************************************************************
 *
 *       Module:    dataio_au4.c
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 26/11/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Reads in the au4 file produced by etpm2au4 which is part of IFP4
 *      Originally written by Paul Eichel, Sandia National Labs
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  19/07/2013
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

#include "sarclib.h"

void read_au4(const char *fname, AuxInfo * aux_info, int num_azi)
{
    FILE * fp;
    char	stmp[80];
    int i;
    
    fp = fopen(fname, "r");
    if (fp == NULL) {
        printf("Can't open file %s\n", fname);
        perror("AU4 open error");
        exit(88);
    }
    
    // skip header down to delimiter string
    //
    for(i=0; i<64; i++)
    {
        fgets(stmp,80,fp);
        if(strncmp(stmp,"**********",10) == 0) break;
    }
    
    for(i = 0; i < num_azi; i++)
    {
        if( fscanf(fp,"%d",&aux_info[i].pulse_number) == EOF) break;
        if( fscanf(fp,"%lf %lf %lf",&aux_info[i].tx_x,
                   &aux_info[i].tx_y,&aux_info[i].tx_z) == EOF) break;
        if( fscanf(fp,"%lf %lf %lf",&aux_info[i].rx_x,
                   &aux_info[i].rx_y,&aux_info[i].rx_z) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].time) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].delta_r) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].fscale) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].c0) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].c1) == EOF) break;
        if( fscanf(fp,"%lf",&aux_info[i].c2) == EOF) break;
    }
    
    fclose(fp);
}

