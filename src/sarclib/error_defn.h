/***************************************************************************
 * 
 *           Module :  error_defn.h
 *          Program :  sarclib
 *       Created by :  Emma Griffiths on 14/12/2004
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *   This file contains the error definitions used in sarclib
 *
 *   DEBUG Defnitions:
 *      Debug >= 10 : minor print output
 *      Debug >= 20 : line and file information on debug output
 *      Debug >= 30 : low level macro debug information - everything
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


#ifndef sarclib_ERRORDEFN_H__
#define sarclib_ERRORDEFN_H__

#include <errno.h>

typedef struct {
    int status;
    int debug;
} SPStatus;

#define NO_ERROR                        (0)
#define NULL_STATUS_POINTER             (1)
#define NULL_POINTER                    (2)
#define NULL_IMAGE                      (3)
#define NULL_STRUCT                     (4)
#define X_PIXEL_SPACING_NOT_EQUAL       (5)
#define Y_PIXEL_SPACING_NOT_EQUAL       (6)
#define NULL_DATA                       (7)
#define INPUT_NX_MISMATCHED             (8)
#define INPUT_NY_MISMATCHED             (9)
#define OUTPUT_NX_MISMATCHED            (10)
#define OUTPUT_NY_MISMATCHED            (11)
#define INPUT_SIZES_MISMATCHED          (12)
#define OUTPUT_SIZES_MISMATCHED         (13)
#define INVALID_FFT_SIZE                (14)
#define NUM_BYTES_WRITTEN_INCORRECT     (15)
#define NULL_FILE_POINTER               (16)
#define ARRAY_SIZES_WRONG               (17)
#define MATRIX_WIDER_THAN_VECTOR        (18)
#define TYPES_NOT_SAME                  (19)
#define INVALID_TYPE                    (20)
#define INVALID_CAST_MODE               (21)
#define LIBRARY_ALREADY_INITIALISED     (22)
#define LIBRARY_NOT_INITIALISED         (23)
#define UNKNOWN_IMAGE                   (24)
#define BAD_IMAGE                       (25)
#define BAD_FILE                        (26)
#define UNSUPPORTED_OPERATOR            (27)
#define ENDIAN_ERROR                    (28)
#define OUT_OF_BOUND                    (29)
#define INTERPOLATION_ERROR             (30)
#define OUT_OF_MEMORY                   (31)
#define CPHD_FILE_INVALID               (32)
#define UNKNOWN_GEOALT_REFERENCE        (33)
#define INVALID_PULSE_INDEX             (34)
#define CPHD_CHANNEL_ERROR              (35)
#define CPHD_DATATYPE_ERROR             (36)
#define CPHD_MAGICNUM_ERROR             (37)
#define CPHD_READ_ERROR                 (38)
#define CPHD_WRITE_ERROR                (39)
#define XML_READ_ERROR                  (40)

/// Check the status flag used to control error reporting in sarclib
/// stat is a pointer to a SPStatus structure
///
#define CHECK_STATUS(stat) {                                                                        \
    if ((stat) == NULL){                                                                            \
        fprintf(stderr, "Null status pointer passed to %s:%d\n", __FILE__, __LINE__);               \
        exit(11);                                                                                   \
    }                                                                                               \
    if ((stat)->status != NO_ERROR){                                                                \
        fprintf(stderr, "Bad status (%d) passed to %s:%d\n",(stat)->status, __FILE__,__LINE__);     \
        exit((stat)->status);                                                                       \
    }                                                                                               \
}

/// Check the status flag used to control error reporting in sarclib. This version is
/// for use when stat is not a pointer
/// stat is a SPStatus structure
///
#define CHECK_STATUS_NON_PTR(stat) {                                                                \
    if (stat.status != NO_ERROR){                                                                   \
        fprintf(stderr, "Bad status (%d) passed to %s:%d\n",stat.status, __FILE__,__LINE__);        \
    }                                                                                               \
}

/// Check to see if a pointer has been initialised
///  point is a pointer
///  s is a pointer to a SPStatus structure
///
#define CHECK_PTR(point,s) {                                                                        \
    if (point == NULL){                                                                             \
        printf("Null pointer found at (%s:%d)\n", __FILE__, __LINE__);                              \
        (s)->status = NULL_POINTER; CHECK_STATUS(s);                                                \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Check an image container to see if it has been correctly initialised
///  imcont is a pointer to a SPImage structure (defined in image.h)
///  s is a pointer to a SPStatus structure
///
#define CHECK_IMCONT(imcont,s) {                                                                    \
    if ((imcont->nx == 0) || (imcont->ny == 0)) {                                                   \
        printf("Null image found (%s:%d)\n", __FILE__, __LINE__);                                   \
        s->status = NULL_IMAGE; return(s);                                                          \
    }                                                                                               \
}

/// Check that two images have the same sample spacing
///
#define CHECK_SPC(img1,img2,s) {                                                                    \
    if ((img1->xspc) != (img2->xspc)){                                                              \
        printf("x-pixel spacing not the same (%s:%d)\n", __FILE__, __LINE__);                       \
        s->status = X_PIXEL_SPACING_NOT_EQUAL; return(s);                                           \
    }                                                                                               \
    if ((img1->yspc) != (img2->yspc)){                                                              \
        printf("y-pixel spacing not the same (%s:%d)\n", __FILE__, __LINE__);                       \
        s->status = Y_PIXEL_SPACING_NOT_EQUAL; return(s);                                           \
    }                                                                                               \
}

/// Check image integrity
/// integ is a pointer to an SPImage. s is a pointer to SPStatus
///
#define CHECK_IMINTEG(integ,s) {                                                                    \
    CHECK_PTR(integ,s);                                                                             \
    CHECK_IMCONT(integ,s);                                                                          \
    CHECK_TYPE_VALID(integ, s);                                                                     \
    if (integ->data.v == NULL) {                                                                    \
        s->status = NULL_DATA;                                                                      \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Check two images have the same number of columns
/// img1 and img2 are pointers to an SPImage struct
/// s is a pointer to an SPStatus struct
///
#define CHECK_COLS(img1,img2,s) {                                                                   \
    if (((img1)->nx) != ((img2)->nx)){                                                              \
        s->status = INPUT_NX_MISMATCHED;                                                            \
        fprintf(stderr, "No cols mismatched at %s:%d nx1 %ld nx2 %ld\n", __FILE__, __LINE__,        \
            (long)(img1)->nx, (long)(img2)->nx);                                                    \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Check two images have the same number of rows
/// img1 and img2 are pointers to an SPImage struct
/// s is a pointer to an SPStatus struct
///
#define CHECK_ROWS(img1,img2,s) {                                                                   \
    if (((img1)->ny) != ((img2)->ny)){                                                              \
        s->status = INPUT_NY_MISMATCHED;                                                            \
        fprintf(stderr, "No rows mismatched at %s:%d ny1 %ld ny2 %ld\n", __FILE__, __LINE__,        \
            (long)(img1)->ny, (long)(img2)->ny);                                                    \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Check two images have the same sizes
/// img1 and img2 are pointers to an SPImage struct
/// s is a pointer to an SPStatus struct
///
#define CHECK_IMSIZEIN(img1,img2,s) {                                                               \
    CHECK_ROWS(img1,img2,s);                                                                        \
    CHECK_COLS(img1,img2,s);                                                                        \
}

/// Check two images have the same sizes and if they don't
/// print out some useful information about them 
/// img1 and img2 are pointers to an SPImage struct
/// s is a pointer to an SPStatus struct
///
#define CHECK_IMSIZEOUT(img1,img2,s) {                                                              \
    if (((img1)->nx) != ((img2)->nx)){                                                              \
        printf("Images have different x sizes (%s:%d)\n", __FILE__, __LINE__);                      \
        s->status = OUTPUT_NX_MISMATCHED;                                                           \
        printf("Error occurred at %s:%d\n", __FILE__, __LINE__);                                    \
        im_info(img1,s);                                                                            \
        im_info(img2,s);                                                                            \
        return(s);                                                                                  \
    }                                                                                               \
    if (((img1)->ny) != ((img2)->ny)){                                                              \
        printf("Images have different y sizes (%s:%d)\n", __FILE__, __LINE__);                      \
        s->status = OUTPUT_NY_MISMATCHED;                                                           \
        printf("Error occurred at %s:%d\n", __FILE__, __LINE__);                                    \
        im_info(img1,s);                                                                            \
        im_info(img2,s);                                                                            \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Macro for checking whether the number of bytes read or written is correct
/// actualNum is the number returned by fread or fwrite. requestedNum is the
/// number of items requested.
/// s is a pointer to SPStatus
///
#define CHECK_BYTES(actualNum,reqestedNum,s) {                                                      \
    if (s->debug >= 30){                                                                            \
        printf("%ld %ld\n", (long)(reqestedNum), (long)(actualNum));                                \
    }                                                                                               \
    if (actualNum != reqestedNum){                                                                  \
        s->status = NUM_BYTES_WRITTEN_INCORRECT;                                                    \
        if (s->status == NUM_BYTES_WRITTEN_INCORRECT){                                              \
            fprintf(stderr, "Number of bytes incorrect (%ld != %ld) at %s:%d\n",                    \
                (long)(actualNum), (long)(reqestedNum), __FILE__,__LINE__);                         \
            if (errno != 0) {                                                                       \
                perror("File I/O error: ");                                                         \
            }                                                                                       \
        }                                                                                           \
    }                                                                                               \
}

/// Macro to validate a file pointer
/// filePointer is the filepointer to be tested
/// s is a pointer to SPStatus
///
#define CHECK_FP(filePointer,s) {                                                                   \
    if (filePointer == NULL){                                                                       \
        fprintf(stderr, "Null file pointer passed to %s:%d\n", __FILE__,__LINE__);                  \
        if (errno != 0) {                                                                           \
            perror("file error: ");                                                                 \
        }                                                                                           \
        s->status = NULL_FILE_POINTER;                                                              \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Macro to check an array sizes against image sizes
/// img1 and img2 are pointers to SPImage
/// x,y,j,k are integers
/// s is a pointer to SPStatus
///
#define CHECK_ARRAY_SIZES(img1,x,y,img2,j,k,s){                                                     \
    if ((x>=img1->nx) || (y>=img1->ny) || (j>=img2->nx) || (k>=img2->ny)){                          \
        s->status = ARRAY_SIZES_WRONG;                                                              \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Macros to check that a matrix has exactly 3 columns (so that it
/// can be multiplied by a vector)
/// mat is of type SPMatrix
/// s is a pointer to SPStatus
///
#define CHECK_MAT_3(mat,s){                                                                         \
    if ((mat.n) != 3){                                                                              \
        s->status = MATRIX_WIDER_THAN_VECTOR;                                                       \
        return(s);                                                                                  \
    }                                                                                               \
}

/// Macros to check that two image types are identical
/// src and dest are of type SPImage
/// s is a pointer to SPStatus
///
#define CHECK_TYPES_EQUAL(src, dest, s) {                                                           \
    if(src->image_type != dest->image_type) {                                                       \
        fprintf(stderr, "Unequal image types (%s != %s) passed to %s:%d\n",                         \
            itype2string(src->image_type), itype2string(dest->image_type), __FILE__,__LINE__);      \
        s->status = TYPES_NOT_SAME; return(s);                                                      \
    }                                                                                               \
}

/// Macro to check that an image has a valid data type
/// a is a pointer to an SPImage struct
/// s is a pointer to SPStatus
///
#define CHECK_TYPE_VALID(a, s)  {                                                                   \
    if(a->image_type == ITYPE_UNKNOWN){                                                             \
        fprintf(stderr, "Invalid type passed to %s:%d\n", __FILE__,__LINE__);                       \
        s->status = INVALID_TYPE; return(s);                                                        \
    }                                                                                               \
}

/// Macros to check that an image type agrees with an expected value.
/// img is a pointer to an SPImage struct
/// type is an SPImageType
/// s is a pointer to SPStatus
///
#define CHECK_TYPE(img, type, s) {                                                                  \
    if(img->image_type != type){                                                                    \
        fprintf(stderr, "Invalid type passed to %s:%d\n", __FILE__,__LINE__);                       \
        s->status = INVALID_TYPE; return(s);                                                        \
    }                                                                                               \
}

/// Macro to print debug information
/// level is an integer
/// status is a pointer to an SPStatus struct
///
#define TRACE(level,status){                                                                             \
    if ((status)->debug >= level){                                                                       \
        printf("%s:%d\n", __FILE__, __LINE__);                                                      \
    }                                                                                               \
}

#endif
