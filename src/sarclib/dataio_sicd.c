/***************************************************************************
 * 
 *           Module :  dataio_sicd.c
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 1/1/2011
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions needed for reading and writing the
 *      SICD (Sensor Independent Complex Data) format.
 *
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/

#include "dataio_sicd.h"

typedef struct {
  char name[14];
  int len;
} NitfTag;

NitfTag  IMGSUBHDR[]  = { { "IM", 2},
			  {"IID1", 10},
			  {"IDATIM", 14},
			  {"TGTID", 17},
			  {"IID2", 80},
			  {"ISCLAS", 1},
			  {"ISCLSY", 2},
			  {"ISCODE", 11},
			  {"ISCTLH", 2},
			  {"ISREL", 20},
			  {"ISDCTP", 2},
			  {"ISDCDT", 8},
			  {"ISDCXM", 4},
			  {"ISDG", 1},
			  {"ISDGDT", 8},
			  {"ISCLTX", 43},
			  {"ISCATP", 1},
			  {"ISCAUT", 40},
			  {"ISCRSN", 1},
			  {"ISSRDT", 8},
			  {"ISCTLN", 15},
			  {"ENCRYP", 1},
			  {"ISORCE", 42},
			  {"NROWS", 8},
			  {"NCOLS", 8},
			  {"PVTYPE", 3},
			  {"IREP", 8},
			  {"ICAT", 8},
			  {"ABPP", 2},
			  {"PJUST", 1},
			  {"ICORDS", 1},
			  {"IGEOLO", 60},
			  {"NICOM", 1},
			  {"ICOMn", 80},
			  {"IC",  2},
			  {"COMRAT", 4},
			  {"NBANDS", 1},
			  {"XBANDS", 5},
			  {"IREPBANDnn", 2},
			  {"ISUBCATnn", 6},
			  {"IFCnn", 1},
			  {"IMFLTnn", 3},
			  {"NLUTSnn", 1},
			  {"NELUTnn", 6},
			  {"LUTDnnm", -1},
			  {"ISYNC", 1},
			  {"IMODE", 1},
			  {"NBPR", 4},
			  {"NBPC", 4},
			  {"NPPBH", 4},
			  {"NPPBV", 4},
			  {"NBPP", 2},
			  {"IDLVL", 3},
			  {"IALVL", 3},
			  {"ILOC", 10},
			  {"IMAG", 4},
			  {"UDIDL", 5},
			  {"UDOFL", 3},
			  {"UDID", -1},
			  {"IXSHDL", 5},
			  {"IXSOFL", 3},
			  {"IXSHD", -1},
			  {"IMDATOFF", 4},
			  {"BMRLNTH", 2},
			  {"TMRLNTH", 2},
			  {"TPXCDLNTH", 2},
			  {"TPXCD",  2},
			  {"BMRnBNDm", 4},
			  {"TMRnBNDm",  4}};

typedef enum {IMGHEADER_IM, IMGHEADER_IID1, IMGHEADER_IDATIM, IMGHEADER_TGTID, IMGHEADER_IID2, IMGHEADER_ISCLAS, IMGHEADER_ISCLSY, IMGHEADER_ISCODE, 
	      IMGHEADER_ISCTLH, IMGHEADER_ISREL, IMGHEADER_ISDCTP, IMGHEADER_ISDCDT, IMGHEADER_ISDCXM, IMGHEADER_ISDG, IMGHEADER_ISDGDT, IMGHEADER_ISCLTX, 
	      IMGHEADER_ISCATP, IMGHEADER_ISCAUT, IMGHEADER_ISCRSN, IMGHEADER_ISSRDT, IMGHEADER_ISCTLN, IMGHEADER_ENCRYP, IMGHEADER_ISORCE, IMGHEADER_NROWS, 
	      IMGHEADER_NCOLS, IMGHEADER_PVTYPE, IMGHEADER_IREP, IMGHEADER_ICAT, IMGHEADER_ABPP, IMGHEADER_PJUST, IMGHEADER_ICORDS, IMGHEADER_IGEOLO, 
	      IMGHEADER_NICOM, IMGHEADER_ICOMn, IMGHEADER_IC, IMGHEADER_COMRAT, IMGHEADER_NBANDS, IMGHEADER_XBANDS, IMGHEADER_IREPBANDnn, IMGHEADER_ISUBCATnn, 
	      IMGHEADER_IFCnn, IMGHEADER_IMFLTnn, IMGHEADER_NLUTSnn, IMGHEADER_NELUTnn, IMGHEADER_LUTDnnm, IMGHEADER_ISYNC, IMGHEADER_IMODE, IMGHEADER_NBPR, 
	      IMGHEADER_NBPC, IMGHEADER_NPPBH, IMGHEADER_NPPBV, IMGHEADER_NBPP, IMGHEADER_IDLVL, IMGHEADER_IALVL, IMGHEADER_ILOC, IMGHEADER_IMAG, 
	      IMGHEADER_UDIDL, IMGHEADER_UDOFL, IMGHEADER_UDID, IMGHEADER_IXSHDL, IMGHEADER_IXSOFL, IMGHEADER_IXSHD, IMGHEADER_IMDATOFF, IMGHEADER_BMRLNTH, 
	      IMGHEADER_TMRLNTH, IMGHEADER_TPXCDLNTH, IMGHEADER_TPXCD, IMGHEADER_BMRnBNDm, IMGHEADER_TMRnBNDm} ImgHeaderEnum;

NitfTag  NITFHEADER[] = { {"FHDR", 9},
			  {"CLEVEL", 2},
			  {"STYPE", 4},
			  {"OSTAID", 10},
			  {"FDT", 14},
			  {"FTITLE", 80},
			  {"FSCLAS", 1},
			  {"FSCLSY", 2},
			  {"FSCODE", 11},
			  {"FSCTLH", 2},
			  {"FSREL", 20},
			  {"FSDCTP", 2},
			  {"FSDCDT", 8},
			  {"FSDCXM", 4},
			  {"FSDG", 1},
			  {"FSDGDT", 8},
			  {"FSCLTX", 43},
			  {"FSCATP", 1},
			  {"FSCAUT", 40},
			  {"FSCRSN", 1},
			  {"FSSRDT", 8},
			  {"FSCTLN", 15},
			  {"FSCOP", 5},
			  {"FSCPYS", 5},
			  {"ENCRYP", 1},
			  {"FBKGC", 3},
			  {"ONAME", 24},
			  {"OPHONE", 18},
			  {"FL", 12},
			  {"HL", 6},
			  {"NUMI", 3},
			  {"LISHnnn", 6},
			  {"LInnn", 10},
			  {"NUMS", 3},
			  {"LSSHnnn", 4},
			  {"LSnnn", 6},
			  {"NUMX", 3},
			  {"NUMT", 3},
			  {"LTSHnnn", 4},
			  {"LTnnn", 5},
			  {"NUMDES", 3},
			  {"LDSHnnn", 4},
			  {"LDnnn", 9},
			  {"NUMRES", 3},
			  {"LRESHnnn", 4},
			  {"LREnnn", 7},
			  {"UDHDL", 5},
			  {"UDHOFL", 3},
			  {"UDHD", -1},
			  {"XHDL", 5},
			  {"XHDLOFL", 3}};

typedef enum {NITFHEADER_FHDR, NITFHEADER_CLEVEL, NITFHEADER_STYPE, NITFHEADER_OSTAID, NITFHEADER_FDT, NITFHEADER_FTITLE, NITFHEADER_FSCLAS, NITFHEADER_FSCLSY, 
	      NITFHEADER_FSCODE, NITFHEADER_FSCTLH, NITFHEADER_FSREL, NITFHEADER_FSDCTP, NITFHEADER_FSDCDT, NITFHEADER_FSDCXM, NITFHEADER_FSDG, NITFHEADER_FSDGDT, 
	      NITFHEADER_FSCLTX, NITFHEADER_FSCATP, NITFHEADER_FSCAUT, NITFHEADER_FSCRSN, NITFHEADER_FSSRDT, NITFHEADER_FSCTLN, NITFHEADER_FSCOP, NITFHEADER_FSCPYS, 
	      NITFHEADER_ENCRYP, NITFHEADER_FBKGC, NITFHEADER_ONAME, NITFHEADER_OPHONE, NITFHEADER_FL, NITFHEADER_HL, NITFHEADER_NUMI, NITFHEADER_LISHnnn, 
	      NITFHEADER_LInnn, NITFHEADER_NUMS, NITFHEADER_LSSHnnn, NITFHEADER_LSnnn, NITFHEADER_NUMX, NITFHEADER_NUMT, NITFHEADER_LTSHnnn, 
	      NITFHEADER_LTnnn, NITFHEADER_NUMDES, NITFHEADER_LDSHnnn, NITFHEADER_LDnnn, NITFHEADER_NUMRES, NITFHEADER_LRESHnnn, NITFHEADER_LREnnn, 
	      NITFHEADER_UDHDL, NITFHEADER_UDHOFL, NITFHEADER_UDHD, NITFHEADER_XHDL, NITFHEADER_XHDLOFL} NITFHeaderEnum;

static SPStatus * load_sicd_meta(SPImage * a, FILE * fp, SICDmetadata * mdata, SPStatus * status);

// Save a complete image to a file
//
SPStatus * im_save_sicd (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status)
{
    FILE *fp;
    int64_t f;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_IMCONT(a,status);
    CHECK_IMINTEG(a,status);
    
    fp = fopen(fname, "w");
    
    CHECK_FP(fp,status);
    
    CHECK_STATUS(status);
    
    im_save_sicd_header(a, mdata, fp, status);
    CHECK_STATUS(status);
    
    f = fwrite(a->data.v,  im_getsizeoftype(a->image_type), a->nx*a->ny, fp);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    fclose(fp);
    
    return(status);
}

// Initialises the line save info structure so that a file can be
// written a few lines at a time
//
SPStatus * im_save_sicd_line_init(SPImage * im, const char * fname, SPImageLineSaveInfo * info, SICDmetadata * mdata, SPStatus * status)
{
    int64_t temp;
    
    CHECK_STATUS(status);
    CHECK_PTR(im,status);
    CHECK_IMCONT(im,status);
    CHECK_IMINTEG(im,status);
    CHECK_PTR(fname,status);
    CHECK_PTR(info,status);
    
    info->fname = strdup(fname);
    info->fp = fopen(fname, "w+");
    CHECK_FP(info->fp,status);
    
    info->nx = im->nx;
    temp = im->ny;
    im->ny = 0;
    im_save_sicd_header(im, mdata, info->fp, status);
    im->ny = temp;
    
    CHECK_STATUS(status);
    
    return(status);
}

// Saves an image to a file that is already partially written. This is useful for writing
// an image to a file a piece at a time. The output file must be already initialised with
// im_save_sicd_line_init().
//
SPStatus * im_save_sicd_line(SPImage * im, SPImageLineSaveInfo * info, SPStatus * status)
{
    int64_t ny;
    int64_t f;
    
    CHECK_STATUS(status);
    CHECK_PTR(im,status);
    CHECK_IMCONT(im,status);
    CHECK_IMINTEG(im,status);
    CHECK_PTR(info,status);
    
    if (im->nx != info->nx)
    {
        fprintf(stderr, "im_save_line: Trying to write out line of different width! (%ld != %ld)\n", (long)im->nx, (long)info->nx);
        status->status = INPUT_NX_MISMATCHED;
        return(status);
    }
    
    fseek(info->fp, sizeof(int)+sizeof(im->nx), SEEK_SET);
    f = fread(&ny, sizeof(ny), 1, info->fp);
    CHECK_BYTES(1, f, status);
    
    ny += im->ny;
    
    fseek(info->fp, sizeof(int)+sizeof(im->nx), SEEK_SET);
    f = fwrite(&ny, sizeof(ny), 1, info->fp);
    CHECK_BYTES(1, f, status);
    
    fseek(info->fp, 0, SEEK_END);
    f = fwrite(im->data.v, im_getsizeoftype(im->image_type), im->nx * im->ny, info->fp);
    CHECK_BYTES(im->nx * im->ny, f, status);
    
    return(status);
}

// Closes down the line save info structure and file
//
SPStatus * im_save_sicd_line_close(SPImageLineSaveInfo * info, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(info,status);
    CHECK_PTR(info->fp, status);
    
    free(info->fname);
    fclose(info->fp);
    return(status);
}


SPStatus * im_save_sicd_header(SPImage *a, SICDmetadata * mdata, FILE * fp, SPStatus * status)
{
    return(status);
}

static SPStatus * load_sicd_meta(SPImage * a, FILE * fp, SICDmetadata * mdata, SPStatus * status)
{
    int i;
    char buff[128];
    
    mdata->numi = 0;
    
    
    // Read in all the NITF tags up to & including NUMI
    //
    for (i = 0; i <= NITFHEADER_NUMI; i++) {
        fread(buff, sizeof(char), NITFHEADER[i].len, fp);
        buff[NITFHEADER[i].len] = 0;
        printf("Reading %s:   %s\n", NITFHEADER[i].name, buff);
        if (strcmp(NITFHEADER[i].name, "NUMI") == 0) {
            mdata->numi = atoi(buff);
        }
    }
    
    for(i = 0; i < mdata->numi; i++) {
        fread(buff, sizeof(char), NITFHEADER[NITFHEADER_LISHnnn].len, fp);
        buff[NITFHEADER[NITFHEADER_LISHnnn].len] = 0;
        printf("    Reading %s:%d   %s\n", NITFHEADER[NITFHEADER_LISHnnn].name, i+1, buff);
        fread(buff, sizeof(char), NITFHEADER[NITFHEADER_LISHnnn].len, fp);
        buff[NITFHEADER[NITFHEADER_LISHnnn].len] = 0;
        printf("    Reading %s:%d   %s\n", NITFHEADER[NITFHEADER_LISHnnn].name, i+1, buff);
        
    }
    
    fread(buff, sizeof(char), NITFHEADER[NITFHEADER_NUMS].len, fp);
    buff[NITFHEADER[NITFHEADER_NUMS].len] = 0;
    printf("    Reading %s:   %s\n", NITFHEADER[NITFHEADER_NUMS].name, buff);
    
    
    
    return(status);
}

// Loads up just the metadata from file
//
SPStatus * im_load_sicd_metadata (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status)  
{
    FILE *fp;
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    im_init(a, status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_sicd_meta(a, fp, mdata, status);
    
    fclose(fp);
    
    return(status);
}

// Loads up a complete image from file
//
SPStatus * im_load_sicd (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status)  
{
    FILE *fp;
    int64_t f;
    int do_swap = 0;
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    
    im_init(a, status);
    
    fp = fopen(fname, "r");
    
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_sicd_meta(a, fp, mdata, status);
    
    im_create (a, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    
    CHECK_IMCONT(a,status);
    
    f = fread(a->data.v, im_getsizeoftype(a->image_type), a->nx*a->ny, fp);
    CHECK_BYTES(a->nx*a->ny, f, status);
    
    swap_data(a, do_swap);
    
    CHECK_IMINTEG(a, status);
    
    fclose(fp);
    
    return(status);
}

// Loads up a subset of an image from a file
//
SPStatus * im_load_sicd_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status) 
{
    FILE *fp;
    int64_t f;
    int64_t y;
    int do_swap = 0;
    SPImage orig;
    SICDmetadata mdata;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    
    fp = fopen(fname, "r");
    CHECK_FP(fp,status);
    CHECK_STATUS(status);
    
    load_sicd_meta(&orig, fp, &mdata, status);
    
    CHECK_IMCONT(a,status);
    if (a->image_type != orig.image_type) {
        fclose(fp);
        fprintf(stderr, "Original image and subimage are different types!");
        status->status = TYPES_NOT_SAME;
        return status;
    }
    
    fseeko(fp, (ox + oy * orig.nx) * im_getsizeoftype(a->image_type), SEEK_CUR);
    
    for(y = 0; y < a->ny; y++) {
        f = fread(&a->data.i8[im_getsizeoftype(a->image_type)*a->nx*y], im_getsizeoftype(a->image_type), a->nx, fp);
        CHECK_BYTES(a->nx, f, status);
        fseeko(fp, (orig.nx - (ox + a->nx) + ox) * im_getsizeoftype(a->image_type), SEEK_CUR);
    }
    
    swap_data(a, do_swap);
    
    CHECK_IMINTEG(a, status);
    
    fclose(fp);
    
    return(status);
}

