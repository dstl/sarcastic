/** @file********************************************************************
 *
 *       Module:    image.h
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 25/10/2004.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *       This file contains functions needed to operate on image structures.
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

#ifndef sarclib_IMAGE_H__
#define sarclib_IMAGE_H__

#ifndef _REENTRANT
#define _REENTRANT
#endif

#include "sarclib.h"

#define im_create(im, type, nx, ny, sx, sy, status)  im_create_function(im, type, nx, ny, sx, sy, status,  (char *)__FILE__, __LINE__)
#define im_destroy(im, status) im_destroy_function(im, status, (char *)__FILE__, __LINE__)

typedef enum {IM_UNKNOWN_ENDIAN, IM_BIG_ENDIAN, IM_LITTLE_ENDIAN} SPEndianType;

typedef enum {DTYPE_COMPLEX = 0x01, DTYPE_SIMPLE = 0x02, DTYPE_FLOAT = 0x04, DTYPE_INTEGER = 0x08} SPDataType;

typedef enum {ITYPE_UNKNOWN,
    ITYPE_POLAR,
    ITYPE_CMPL_FLOAT,
    ITYPE_FLOAT,  ITYPE_DOUBLE,
    ITYPE_INT64,  ITYPE_INT32,  ITYPE_INT16,  ITYPE_INT8,
    ITYPE_UINT64, ITYPE_UINT32, ITYPE_UINT16, ITYPE_UINT8,
    ITYPE_VECTOR,
    ITYPE_CMPL_INT64, ITYPE_CMPL_INT32, ITYPE_CMPL_INT16, ITYPE_CMPL_INT8, ITYPE_CMPL_NIBBLE } SPImageType;

typedef enum {CAST_MODULUS, CAST_PHASE, CAST_UNIT, SET_MODULUS, SET_PHASE, SET_REAL, SET_IMAG, CAST_REAL, CAST_IMAG} SPCastMode;

typedef enum {X_DIR, Y_DIR, Z_DIR} SPDirection;

typedef struct {
    int64_t nx;                         ///< no columns of data
    int64_t ny;                         ///< no rows of data
    double xspc;                        ///< pixel spacing in x
    double yspc;                        ///< pixel spacing in y
    SPImageType image_type;             ///< What type of image data do it contain
    uint64_t  id;                       ///< Private ID used for image tracking
    void * user_data;                   ///< Generic pointer to let user do interesting stuff!
    union {
        SPCmplxPol   *pol;              ///< pointer to the complex data array - float polar style
        SPCmplx      *cmpl_f;           ///< pointer to the complex data array - float cartesian style
        SPCmplxInt64 *cmpl_i64;         ///< pointer to the complex data array - int64 cartesian style
        SPCmplxInt32 *cmpl_i32;         ///< pointer to the complex data array - int32 cartesian style
        SPCmplxInt16 *cmpl_i16;         ///< pointer to the complex data array - int16 cartesian style
        SPCmplxInt8  *cmpl_i8;          ///< pointer to the complex data array - int8 cartesian style
        SPCmplxNibb  *cmpl_nibb;        ///< pointer to the complex data array - nibble cartesian style
        float        *f;                ///< pointer to float data
        double       *d;                ///< pointer to double data
        int64_t      *i64;              ///< pointer to int64 data
        int32_t      *i32;              ///< pointer to int32 data
        int16_t      *i16;              ///< pointer to int16 data
        int8_t       *i8;               ///< pointer to int8 data
        uint64_t     *ui64;             ///< pointer to uint64 data
        uint32_t     *ui32;             ///< pointer to uint32 data
        uint16_t     *ui16;             ///< pointer to uint16 data
        uint8_t      *ui8;              ///< pointer to uint8 data
        SPVector     *vect;             ///< Pointer to vector data
        void         *v;                ///< Generic pointer to allow artful hacking
    } data;
} SPImage;


#define im_init_status(SPStatVar, debug_lvl) {                  \
    SPStatVar.status = NO_ERROR; SPStatVar.debug = debug_lvl;   \
}

#ifdef __CUFILE__
extern "C"{
#endif
    
    // Functions in image.c
    //
    /// create an image structure 
    SPStatus* im_create_function(SPImage *a, SPImageType t, int64_t sizex, int64_t sizey, double spx, double spy, SPStatus *status, const char * fname, int line);
    
    /// destroy an image structure to free up the memory
    ///
    SPStatus* im_destroy_function (SPImage *a, SPStatus *status, const char *fname, int line);
    
    /// initialise an image structure
    ///
    SPStatus* im_init (SPImage *a, SPStatus *status);
    
    /// clones an image (i.e. copies only the data pointer)
    ///
    SPStatus* im_clone (SPImage *a, SPImage *b, SPStatus *status);
    
    /// creates a copy of an image
    ///
    SPStatus* im_copy (SPImage *orig, SPImage *copy, SPStatus *status);
    
    /// fills array with data, for use inside an x-y loop
    ///
    SPStatus* im_fill (SPImage *a, int64_t x, int64_t y, SPImage *b, int64_t j, int64_t k, SPStatus *status);
    
    /// Function that must must be called before you start to use these functions
    ///
    SPStatus* im_init_lib(SPStatus * status, const char * prog_name, int argc, char * const * argv);
    
    /// Prints out the current status of the im library
    ///
    SPStatus* im_status(SPStatus * status);
    
    /// Closes down the library when you've finished with it
    ///
    SPStatus* im_close_lib(SPStatus * status);
    
    /// This function returns the size of an element of the given ImageType
    ///
    size_t im_getsizeoftype(SPImageType t);
    
    /// prints the line which an image was allocated on
    ///
    void im_info(SPImage *a, SPStatus *status);
    
    /// Returns the endianness of the machine the code is being executed on
    ///
    SPEndianType im_machine_type(void);
    
    /// zeroises the data part of an array
    ///
    SPStatus *im_zero_data(SPImage *in, SPStatus *status);  
    
    /// returns the pixel value as a double - for speed issues,
    /// don't use this in a loop!
    ///
    double im_get_pixel_as_double(SPImage * im, int64_t x, int64_t y, SPStatus * status);
    
    
    // Functions in image_average.c
    //
    
    /// Make a smaller image by averaging a larger one!
    ///
    SPStatus* im_average(SPImage * src, SPImage * dest, SPStatus * status);  
    
    /// Functions in image_convertion.c
    ///
    
    /// cart to polar
    ///
    SPStatus* im_cart_to_polar (SPImage *a, SPImage *out, SPStatus *status);
    
    /// Convert polar data to cartesian
    ///
    SPStatus* im_polar_to_cart(SPImage *a, SPImage * out, SPStatus *status);
    
    /// Convert the type of an image
    ///
    SPStatus* im_cast(SPImage *src, SPCastMode mode, SPImage *dst, SPStatus * status);
    
    /// Converts two images (re & im) to a complex float
    ///
    SPStatus *im_cmplx_f_make(SPImage *real, SPImage *imaginary, SPImage *out, SPStatus *status);  
    
    /// Functions in image_extract.c
    ///
    
    /// extract a bit of an image from larger structure
    ///
    SPStatus* im_extract (SPImage *src, int64_t startx, int64_t starty, SPImage * out, SPStatus *status);
    
    /// shifts the elements of the image circularly
    ///
    SPStatus* im_circshift (SPImage *a, int64_t shift_x, int64_t shift_y, SPStatus *status);
    
    /// Insert the src image into dest image at (x,y)
    ///
    SPStatus* im_insert(SPImage * src, int64_t x, int64_t y, SPImage * dest, SPStatus * status);   

    
    /// Functions in image_fft.c
    ///
    
    /// fft the contents of an image
    ///
    SPStatus* im_fft(SPImage *a, FFTMODE mode, SPStatus *status);  
    
    /// Functions in image_fftw.c
    ///
    
    /// FFT of a real 1D to complex 1D image
    SPStatus* im_fftw_r2c(SPImage * inp, FFTMODE mode, SPImage * out, SPStatus * status);

    /// fft the contents of an image
    ///
    SPStatus* im_fftw(SPImage *a, FFTMODE mode, SPStatus *status);  
    
    /// Functions in image_math.c
    ///
    
    /// adds the data part of two image structures
    ///
    SPStatus* im_add (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// subtracts the data part of two image structures
    ///
    SPStatus* im_sub (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// multiplies the data part of two image structures
    ///
    SPStatus* im_mult (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// divides the data part of two image structures
    ///
    SPStatus* im_div (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// complex conjugate of the data part of an image structure
    ///
    SPStatus* im_conjg (SPImage *a, SPImage *out, SPStatus *status);
    
    /// multiplies the data part of an image structure by a scalar
    ///
    SPStatus* im_mult_scalar(SPImage *im, SPImageType c_type, void *con, SPStatus *status);
    
    /// as above but outputs image
    ///
    SPStatus* im_mult_scalar_o(SPImage *im, SPImageType c_type, void *con, SPImage *out, SPStatus *status);
    
    /// subtracts each row of a col of doubles to each row of the data part of an image structure
    ///
    SPStatus* im_colsub (SPImage *a, double *col, SPImage *out, SPStatus *status);
    
    /// adds each row of a column of doubles to each row of the data part of an image structure
    ///
    SPStatus* im_coladd (SPImage *a, double *col, SPImage *out, SPStatus *status);
    
    /// correlate the data part of two image structures
    ///
    SPStatus* im_corr (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// add up each column in the image
    ///
    SPStatus* im_sumcols(SPImage *a, SPImage *out, SPStatus *status);
    
    /// add up each row in the image
    ///
    SPStatus* im_sumrows(SPImage *a, SPImage *out, SPStatus *status);
    
    /// adds a constant onto the supplied image
    ///
    SPStatus* im_add_scalar(SPImage * im, SPImageType c_type, void * con, SPStatus * status);
    
    /// as above and outputs image
    ///
    SPStatus* im_add_scalar_o(SPImage * im, SPImageType c_type, void * con, SPImage *out, SPStatus * status);
    
    /// Returns the next highest power of two
    ///
    int64_t next_power_2(int64_t num);
    
    /// Multiplies (pointwise) a 2D image by a 1D image.
    ///
    SPStatus* im_mult_line (SPImage *a, SPImage *b, SPDirection dir, SPImage *out, SPStatus *status);
    
    /// Gives the dot product of each element of two vector images, outputs it as a float
    ///
    SPStatus* im_dot (SPImage *a, SPImage *b, SPImage *out, SPStatus *status);
    
    /// Gives the fmod (# man fmod) of each element of an image and outputs it as a double.
    ///
    SPStatus* im_fmod (SPImage *a, double div, SPImage *out, SPStatus *status);
    
    /// Gives the sin (# man sin) of each element of an image and outputs it as a double.
    ///
    SPStatus* im_sin (SPImage *in, SPImage *out, SPStatus *status);
    
    /// Gives the cos (# man cos) of each element of an image and outputs it as a double.
    ///
    SPStatus* im_cos (SPImage *in, SPImage *out, SPStatus *status);
    
    /// Gives the rint (# man rint) of each element of an image and outputs it as a int64_t.
    ///
    SPStatus* im_rint (SPImage *a, SPImage *out, SPStatus *status);
    
    // Functions in image_stats.c
    //
    /// Computes some stats about the image
    ///
    SPStatus* im_stats(SPImage * a, SPStatsInfo * info, SPStatus * status);
    
    /// Create a histogram of the image
    ///
    SPStatus* im_histo(SPImage * a, double min_bin, double max_bin, HistoMode mode,
                       int64_t num_bins, SPImage * hist, SPStatus * status);  
    
    /// Functions in image_transpose.c
    ///
    /// transpose of an image structure
    ///
    SPStatus* im_transp (SPImage *a, SPStatus *status);  
    
    
    // Functions in image_pad.c
    //
    /// zero pads an image, placing zeros on right and bottom of image
    ///
    SPStatus* im_zpad(SPImage *a, int64_t new_size_x, int64_t new_size_y, SPImage *out, SPStatus *status);  
    
    // Functions in image_recip. c
    //
    /// This is a function to take the reciprocal of an image, in and out can be the same image if necessary
    ///
    SPStatus * im_recip(SPImage *in, SPImage *out, SPStatus *status);
    
    // Functions in image_unwrap.c
    //
    /// unwrap phase from the supplied 1D array
    ///
    SPStatus * im_unwrap (SPImage *phas, SPImage * unwrapped, SPStatus *status);  
    
    /// Perform matrix multiplication
    ///
    SPStatus * im_matrix_mult(SPImage *in1, SPImage *in2, SPImage *out, SPStatus * status);
    
    // Functions in image_random.c
    //
    /// Fill an array with gauss noise
    ///
    SPStatus * im_gauss(SPImage * inp, void * mean, void * sd, SPStatus * status);
    
    /// Set the seed of the PRNG
    ///
    SPStatus * im_init_rand(long seed, SPStatus * status);
    
    /// Generate a random number - with gaussian distrubution
    ///
    double rand_gauss(void); 
    
#ifdef __CUFILE__
}
#endif


#endif
