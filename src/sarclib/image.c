/***************************************************************************
 * 
 *           Module :  image.c
 *          Program :  sarclib
 *       Created by :  Emma Griffiths on 25/10/2004
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *       This file contains functions needed to operate on image structures.
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

#include <unistd.h>
#include "sarclib.h"

#include <fftw3.h>

#define MAX_FFTW_WIDSOM_FNAME_LEN  (256)

// This defines a structure which is used as a linked list to store what images
// have been allocated so we can easily check for memory leaks. It won't
// take much to extend im_calloc & im_destroy to do much more such as
// integrety checking
//
struct memList_ {
    struct memList_ * next;
    void *   data;
    size_t   nelm;
    size_t   elm_size;
    int64_t  nx;
    int64_t  ny;
    char  *  fname;
    int      line;
    int      active;
    uint64_t id;
};

typedef struct memList_ memList;

// This is the linked list of information about the images which we allocated
//
static memList * all_mem = NULL;

// This is just a count so we can figure out what image we are trying to destroy if its
// been cloned or copied etc..
//
static uint64_t count = 1;

// This is a mutex lock so that we don't try and read and write to the above linked list
// and so we increment count in a controlled manner. This is only necessary in a pthread'ed
// program, but we do it just in case.
//
pthread_mutex_t memlist_lock;

// This is a mutex lock for the fftw_plan function
//
pthread_mutex_t fftw_plan_lock;

// This contains a copy of the program name
//
char * program_name = NULL;

// These integers come from the getopt function
//
extern int optind, opterr, optopt;

// This will contain the input parameter file name to be used by the input_* functions
//
char * input_parameter_fname = NULL;

// This contains the output parameter file name to be used by the input_* functions
//
char * output_parameter_fname = NULL;

// Should we ask the user questions - default is we ask questions
//
int no_ask_questions = FALSE;

// Should we rename the parameter file at the close of the program to be $PWD/.<program_name>rc
//
int no_parameter_file_rename = FALSE;

// This is the prototype for the private function im_calloc - more comments below
// It allocates mem for image and maintains linked list of image information
//
static void * im_calloc(int64_t nx, int64_t ny, size_t elm_size, uint64_t * id, const char * fname, int line, SPStatus * status);

// This is the low-level function that allocates all the memory for the supplied image.
// This is called from the im_create macro (defined in image.h) with just __FILE__ &
// __LINE__ tacked on so that the linked list above knows from where the im_create
// was actually called.
//
// Need to supply the image type, the size in x & y, the spacing in x & y and as ever,
// the status.
//
SPStatus* im_create_function (SPImage *a, SPImageType t, int64_t sizex, int64_t sizey, double spx, double spy, SPStatus *status, const char *fname, int line)
{
    if (status->debug >= 20)
    {
        printf("%s:%d - %s:%d\n", fname, line, __FILE__, __LINE__);
        printf("%p\n", a);
    }
    CHECK_STATUS(status);
    if (status->debug >= 20)
    {
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("%p\n", a);
    }
    im_init(a, status);
    
    a->nx = sizex;
    a->ny = sizey;
    a->xspc = spx;
    a->yspc = spy;
    a->image_type = t;
    
    CHECK_IMCONT (a,status);
    
    a->data.v = im_calloc(a->nx, a->ny, im_getsizeoftype(t), &a->id, fname, line, status);
    
    CHECK_PTR(a->data.v, status);
    CHECK_STATUS(status);
    
    return(status);
}

// This function destroys (free's up memory assoc'ed with an image. Like the
// create function above, its called from a helper macro that tacks on the __FILE__ and
// __LINE__.
//
// The function goes through the linked list of image information, updates that information then
// free's the memory. All while handling the mutex lock as required.
//
SPStatus* im_destroy_function (SPImage *a, SPStatus *status, const char *fname, int line)
{
    memList * i;
    
    CHECK_STATUS (status);
    
    // Check to see if we are destroying an image which has only been init'ed
    //
    if (a->nx == 0 && a->ny == 0 && a->id == 0 && a->data.v == NULL && a->image_type == ITYPE_UNKNOWN)
    {
        return status;
    }
    
    pthread_mutex_lock(&memlist_lock);
    
    for (i = all_mem; i->next != NULL && i->id != a->id; i = i->next);
    if (i->id == 0)
    {
        fprintf(stderr, "Trying to destroy an image with an ID (%lx) which is not in the list!  @ %s:%d\n", (long unsigned int)a->id, fname, line);
        status->status = UNKNOWN_IMAGE;
        pthread_mutex_unlock(&memlist_lock);
        return(status);
    }
    
    if (i->data != a->data.v && i->active == 1)
    {
        fprintf(stderr, "Trying to destroy an image (%16lx) with a data ptr that has changed! (%p != %p) @ %s:%d\n", (long unsigned int)a->id, i->data, a->data.v, fname, line);
        status->status = BAD_IMAGE;
        pthread_mutex_unlock(&memlist_lock);
        return(status);
    }
    
    if (i->active == 0)
    {
        fprintf(stderr, "Trying to free an image which has already been flagged as inactive! @ %s:%d\n", fname, line);
        status->status = BAD_IMAGE;
        pthread_mutex_unlock(&memlist_lock);
        return(status);
    }
    
    i->active = 0;
    i->id = 0;
    free (a->data.v);
    i->data = NULL;
    free(i->fname);
    i->fname = NULL;
    
    pthread_mutex_unlock(&memlist_lock);
    
    im_init(a, status);
    return(status);
}

// Sets the members of the image structure to some sensible values
//
SPStatus* im_init (SPImage *a, SPStatus *status)  
{
    CHECK_STATUS(status);
    if (status->debug >= 20)
    {
        printf("%s:%d : Ptr to image being initialised: %p\n", __FILE__, __LINE__,a);
    }
    
    a->nx = 0;
    a->ny = 0;
    a->data.v = NULL;
    a->xspc = 1.0;
    a->yspc = 1.0;
    a->image_type = ITYPE_UNKNOWN;
    a->id = 0;
    
    TRACE(25, status);
    status->status = NO_ERROR;
    return(status);
}


// Create an exact copy of an image - note that you have two images which *share* the same
// data, so if you destroy one, you can't get at the data using the other image.
//
SPStatus* im_clone (SPImage *a, SPImage *b, SPStatus *status) 
{
    TRACE(25, status);
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_PTR(b,status);
    CHECK_IMCONT(a,status);
    CHECK_IMINTEG(a,status);
    
    b->nx = a->nx;
    b->ny = a->ny;
    b->xspc = a->xspc;
    b->yspc = a->yspc;
    b->data.v = a->data.v;
    b->image_type = a->image_type;
    b->id = a->id;
    TRACE(25, status);
    return(status);
}

// Makes an identical copy of the original supplied image
//
SPStatus* im_copy (SPImage *orig, SPImage *copy, SPStatus *status)
{
    TRACE(25, status);
    CHECK_STATUS(status);
    CHECK_PTR(orig,status);
    CHECK_IMCONT(orig,status);
    CHECK_IMINTEG(orig,status);
    
    if (copy->data.v == NULL && copy->nx == 0 && copy->ny == 0)
    {
        im_create(copy, orig->image_type, orig->nx, orig->ny, orig->xspc, orig->yspc, status);
    }
    
    CHECK_IMSIZEOUT(orig, copy, status);
    
    CHECK_STATUS(status);
    
    if (copy->nx == orig->nx && copy->ny == orig->ny && copy->image_type == orig->image_type)
    {
        memcpy(copy->data.v, orig->data.v, im_getsizeoftype(orig->image_type) * orig->nx * orig->ny);
    }
    else
    {
        status->status = ARRAY_SIZES_WRONG;
    }
    TRACE(25, status);
    return(status);
}


// Copies a pixel from b(j,k) to a(x,y)
//
SPStatus* im_fill (SPImage *a, int64_t x, int64_t y, SPImage *b, int64_t j, int64_t k, SPStatus *status)
{
    CHECK_ARRAY_SIZES(a,x,y,b,j,k,status);
    CHECK_STATUS(status);
    CHECK_TYPES_EQUAL(a, b, status);
    TRACE(25, status);
    
    switch (a->image_type)
    {
        case ITYPE_POLAR:
            a->data.pol[x+a->nx*y] = b->data.pol[j+b->nx*k];
            break;
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_fill\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    TRACE(25, status);
    return(status);
}

// Allocates the memory for a new image. Also maintains the linked list of image information.
// This is a PRIVATE function and shouldn't be called by anyone other than the library.
// The image size (nx & ny), the size of the element, the unique ID, and the filename and
// line number.
//
static void * im_calloc(int64_t nx, int64_t ny, size_t elm_size, uint64_t *id, const char * fname, int line, SPStatus * status)
{
    memList * new = NULL;
    memList * i;
    
    pthread_mutex_lock(&memlist_lock);
    
    for (i = all_mem; i != NULL && new == NULL; i = i->next)
    {
        if (i->active == 0)
        {
            new = i;
        }
    }
    
    if (new == NULL)
    {
        new = calloc(1, sizeof(memList));
        
        if (new == NULL)
        {
            status->status = NULL_POINTER;
            fprintf(stderr, "Failed to calloc mem for new memory list element!\n");
            return(NULL);
        }
        new->next = all_mem;
        all_mem = new;
    }
    
    new->nelm = nx * ny;
    new->elm_size = elm_size;
    new->nx = nx;
    new->ny = ny;
    new->fname = strdup(fname);
    new->line = line;
    new->active = 0;
    new->data = calloc(nx*ny, elm_size);
    
    if (new->data == NULL)
    {
        status->status = NULL_POINTER;
	fprintf(stderr, "Failed to calloc mem for new image data (%"PRId64" x %"PRId64" of %ld bytes requested)!\n", nx, ny, elm_size);
        return(NULL);
    }
    new->active = 1;
    
    new->id = count;
    *id = count;
    count++;
    
    pthread_mutex_unlock(&memlist_lock);
    
    return (new->data);
}

// This function initialises the library and sets status to be OK.
//
// It sets up the linked list of image information, the mutex lock and count.
// It copies the program name and can parse the cmd line options.
//
// This function is not multithread safe, and must be called once, and only
// once, in a program.
//
SPStatus* im_init_lib(SPStatus * status, const char * progname, int argc, char * const *argv)
{
    char * home = getenv("HOME");
    size_t len;
    FILE * fp;
    int c;
    int no_read_input;
    
    status->status = NO_ERROR;
    CHECK_STATUS(status);
    if (all_mem != NULL)
    {
        status->status = LIBRARY_ALREADY_INITIALISED;
        fprintf(stderr, "Library already initalised!\n");
        return (status);
    }
    
    all_mem = calloc(1, sizeof(memList));
    
    if (all_mem == NULL)
    {
        status->status = NULL_POINTER;
        fprintf(stderr, "Failed to calloc mem for memory list element!\n");
        return(status);
    }
    pthread_mutex_init(&memlist_lock, NULL);
    count = 1;
    
    all_mem->next = NULL;
    all_mem->active = 0;
    
    program_name = strdup(progname);
    if (!program_name)
    {
        status->status = NULL_POINTER;
        fprintf(stderr, "Failed to calloc mem for program_name!\n");
        return(status);
    }
    
    no_read_input = FALSE;
    
    do {
        c = getopt(argc, argv, "s:r:dh?n");
        if (c != -1)
        {
            switch(c)
            {
                case 's':
                    output_parameter_fname = strdup(optarg);
                    if (output_parameter_fname == NULL)
                    {
                        fprintf(stderr, "Failed to strdup output_parameter_fname!\n");
                        exit(8881);
                    }
                    no_parameter_file_rename = TRUE;
                    
                    break;
                case 'r':
                    input_parameter_fname = strdup(optarg);
                    if (input_parameter_fname == NULL)
                    {
                        fprintf(stderr, "Failed to strdup input_parameter_fname!\n");
                        exit(8882);
                    }
                    break;
                case 'd':
                    no_read_input = TRUE;
                    break;
                case 'n':
                    no_ask_questions = TRUE;
                    break;
                case 'R':
                    no_parameter_file_rename = TRUE;
                    break;
                case 'h':
                case '?':
                    printf("Usage:\n");
                    printf("-s <filename>   save input parameters to file <filename>, assumes -R as well\n");
                    printf("-r <filename>   read input parameters from file <filename>\n");
                    printf("-d              don't take any defaults from a file - use the ones supplied by the program\n");
                    printf("-R              don't rename the output parameter file to $PWD/.%src\n", program_name);
                    printf("-n              don't ask questions, just use the defaults\n");
                    printf("-h              print out this help\n");
                    printf("-?              print out this help\n");
                    printf("\n\nsarclib version %s\n", FULL_VERSION);
                    break;
            }
        }
    } while (c != -1);
    
    if (home != NULL)
    {
        len = strlen(progname) + strlen(home) + 40;
    }
    else
    {
        len = strlen(progname) + 40;
    }
    
    if (!no_read_input)
    {
        if (input_parameter_fname == NULL)
        {
            input_parameter_fname = calloc(len, sizeof(char));
            snprintf(input_parameter_fname, len, ".%src",progname);
            fp = fopen(input_parameter_fname, "r");
            
            if (fp == NULL && home != NULL)
            {
                snprintf(input_parameter_fname, len, "%s/.%src", home, progname);
                fp = fopen(input_parameter_fname, "r");
            }
            
            if (fp == NULL)
            {
                free(input_parameter_fname);
                input_parameter_fname = NULL;
                printf("Cannot read parameter file from either CWD or $HOME - using program defaults\n");
            }
            else
            {
                fclose(fp);
            }
        }
        else
        {
            // Because this filename has been requested by the user, shout and then die if we can't open that file
            //
            fp = fopen(input_parameter_fname, "r");
            if (fp == NULL)
            {
                printf("Failed to open file \"%s\" for read\n", input_parameter_fname);
                exit(3489);
            }
            fclose(fp);
        }
    }
    
    // We test open for write (rather than append like we actually open the file for the input_* routines), as that
    // has the benefit of emptying the file
    //
    if (output_parameter_fname == NULL)
    {
        output_parameter_fname = calloc(len, sizeof(char));
        snprintf(output_parameter_fname, len, ".%src_%d",progname, (int)getpid());
        fp = fopen(output_parameter_fname, "w");
        if (fp == NULL && home != NULL)
        {
            snprintf(output_parameter_fname, len, "%s/.%src_%d", home, progname, (int)getpid());
            fp = fopen(output_parameter_fname, "w");
        }
        
        if (fp == NULL)
        {
            free(output_parameter_fname);
            output_parameter_fname = NULL;
        }
        else
        {
            fclose(fp);
        }
    }
    else
    {
        fp = fopen(output_parameter_fname, "w");
        if (fp == NULL)
        {
            printf("Failed to open file \"%s\" for write\n", output_parameter_fname);
            exit(3489);
        }
        fclose(fp);
    }
    
    char fftw_wisdom_file[MAX_FFTW_WIDSOM_FNAME_LEN];

    char * hm = getenv("HOME");
    if (hm == NULL) {
      snprintf(fftw_wisdom_file, MAX_FFTW_WIDSOM_FNAME_LEN, ".fftw_wisdom");
    } else {
      snprintf(fftw_wisdom_file, MAX_FFTW_WIDSOM_FNAME_LEN, "%s/.fftw_wisdom", hm);
    }
    fftwf_import_wisdom_from_filename(fftw_wisdom_file);

    pthread_mutex_init(&fftw_plan_lock, NULL);
    
    return (status);
}

// This function returns a list of all the currently allocated images.
//
SPStatus* im_status(SPStatus * status)
{
    memList * i;
    CHECK_STATUS(status);
    int64_t count = 0;
    int tcount = 0;
    
    if (all_mem == NULL)
    {
        status->status = LIBRARY_NOT_INITIALISED;
        fprintf(stderr, "Library does not appear to be initalised. Did you call im_init_lib()?\n");
        return(status);
    }
    
    pthread_mutex_lock(&memlist_lock);
    for (i = all_mem; i != NULL; i = i->next) {
        if (i->active == 1) {
            tcount++;
        }
    }
    
    // Alter it so that it only prints out stuff if we have things allocated & active or if you have a
    // non-zero debug level
    //
    if (tcount > 0 || status->debug > 0) {
        printf("The following images are still allocated:\n");
        for (i = all_mem; i != NULL; i = i->next) {
            if (i->active == 1) {
                printf("File:%20s  Line:%5d  Nx:%8ld  Ny:%8ld  ID: %016lx   ptr:%p\n", i->fname, i->line, (long)i->nx, (long)i->ny, (long unsigned int)i->id, i->data);
            } else {
                count++;
            }
        }
        printf("\nThere are %ld unused nodes in the memory linked list\n", (long)count);
    }
    
    pthread_mutex_unlock(&memlist_lock);
    
    return (status);
}

// Shuts down the library.
// At moment just prints out a list of active imagess - probably should destroy them.
//
SPStatus* im_close_lib(SPStatus * status)
{
    char * tmp;
    im_status(status);
    
    if (!no_parameter_file_rename) {
        tmp = calloc(strlen(program_name)+20, sizeof(char));
        if (tmp == NULL)
        {
            fprintf(stderr, "Failed to calloc array for rename in close_lib!!\n");
            exit(6655);
        }
        snprintf(tmp, strlen(program_name)+20, ".%src", program_name);
        if (rename(output_parameter_fname, tmp) == -1)
        {
            perror("Rename failed");
        }
        free(tmp);
    }

    char fftw_wisdom_file[MAX_FFTW_WIDSOM_FNAME_LEN];

    char * hm = getenv("HOME");
    if (hm == NULL) {
      snprintf(fftw_wisdom_file, MAX_FFTW_WIDSOM_FNAME_LEN, ".fftw_wisdom");
    } else {
      snprintf(fftw_wisdom_file, MAX_FFTW_WIDSOM_FNAME_LEN, "%s/.fftw_wisdom", hm);
    }
    fftwf_export_wisdom_to_filename(fftw_wisdom_file);
    
    return (status);
}

// This function returns the size of an element of the given ImageType.
//
size_t im_getsizeoftype(SPImageType t)
{
    size_t sz;
    
    switch(t)
    {
        case ITYPE_POLAR:
            sz = sizeof(SPCmplxPol);
            break;
            
        case ITYPE_CMPL_FLOAT:
            sz = sizeof(SPCmplx);
            break;
            
        case ITYPE_FLOAT:
            sz = sizeof(float);
            break;
            
        case ITYPE_DOUBLE:
            sz = sizeof(double);
            break;
            
        case ITYPE_INT64:
            sz = sizeof(int64_t);
            break;
            
        case ITYPE_INT32:
            sz = sizeof(int32_t);
            break;
            
        case ITYPE_INT16:
            sz = sizeof(int16_t);
            break;
            
        case ITYPE_INT8:
            sz = sizeof(int8_t);
            break;
            
        case ITYPE_UINT64:
            sz = sizeof(uint64_t);
            break;
            
        case ITYPE_UINT32:
            sz = sizeof(uint32_t);
            break;
            
        case ITYPE_UINT16:
            sz = sizeof(uint16_t);
            break;
            
        case ITYPE_UINT8:
            sz = sizeof(uint8_t);
            break;
            
        case ITYPE_VECTOR:
            sz = sizeof(SPVector);
            break;
            
        case ITYPE_CMPL_INT64:
            sz = sizeof(SPCmplxInt64);
            break;
            
        case ITYPE_CMPL_INT32:
            sz = sizeof(SPCmplxInt32);
            break;
            
        case ITYPE_CMPL_INT16:
            sz = sizeof(SPCmplxInt16);
            break;
            
        case ITYPE_CMPL_INT8:
            sz = sizeof(SPCmplxInt8);
            break;
            
        default:
            fprintf(stderr, "Unknown/Invalid type passed to im_getsizeoftype (%d)\n", t);
            sz = 0;
            break;
    }
    
    return (sz);
}

void im_info(SPImage *a, SPStatus *status)
{
    
    memList *i;
    
    pthread_mutex_lock(&memlist_lock);
    
    for (i = all_mem; i->next != NULL && i->id != a->id; i = i->next);
    
    if (i->id == 0)
    {
        fprintf(stderr, "Unknown image ID\n");
        pthread_mutex_unlock(&memlist_lock);
        return;
    }
    
    if (i->data != a->data.v && i->active == 1)
    {
        fprintf(stderr, "Image ID has a data ptr that has changed! (%p != %p) @ %s:%d\n", i->data, a->data.v, i->fname, i->line);
        pthread_mutex_unlock(&memlist_lock);
        return;
    }
    
    if (i->active == 0)
    {
        fprintf(stderr, "This image is not active! @ %s:%d\n", i->fname, i->line);
        pthread_mutex_unlock(&memlist_lock);
        return;
    }
    
    printf("This image was allocated at %s : %d\n", i->fname, i->line);
    printf("The size of the image is : nx %ld ny %ld\n", (long)a->nx, (long)a->ny);
    
    pthread_mutex_unlock(&memlist_lock);
}

// Returns the endianness of the machine the code is being executed on
//
SPEndianType im_machine_type(void)  
{
    SPEndianType ret = IM_UNKNOWN_ENDIAN;
    union {char c[4]; int i;} f;
    
    // This lovely little hack makes a union of 4 chars and an int. It sets the int equal to one, and
    // then just sees which of the chars (either the first or the last) ends up being set to one - which
    // one it is depends on if its big or little endian!!
    //
    f.i = 1;
    
    if (f.c[0] == 1) {
        ret = IM_LITTLE_ENDIAN;
    }
    
    if (f.c[3] == 1) {
        ret = IM_BIG_ENDIAN;
    }
    
    return ret;
}

// zeroises the data part of an array
//
SPStatus *im_zero_data(SPImage *in, SPStatus *status) 
{
    CHECK_IMINTEG(in, status);
    CHECK_STATUS(status);
    
    memset(in->data.v, 0, im_getsizeoftype(in->image_type)*in->nx*in->ny);
    
    return(status);
}


#define GET_PIX(t, x, y) ((double) im->data. t [x + im->nx * y])

double im_get_pixel_as_double(SPImage * im, int64_t x, int64_t y, SPStatus * status)
{
    double val = -1;
    
    if (status->status != NO_ERROR) {
        fprintf(stderr, "Bad status passed to im_get_pixel_as_double.\n");
        return val;
    }
    
    if ((im == NULL) || (im->nx == 0) || (im->ny == 0)) {
        fprintf(stderr, "Null image passed to im_get_pixel_as_double.\n");
        status->status = NULL_IMAGE;
        return val;
    }
    
    if (x >= im->nx || x < 0) {
        fprintf(stderr, "Trying to access x pixel out of range (%"PRId64"  is <0 or >=%"PRId64") in im_get_pixel_as_double\n", x, im->nx);
        status->status = OUT_OF_BOUND;
        return val;
    }
    
    if (y >= im->ny || y < 0) {
        fprintf(stderr, "Trying to access y pixel out of range (%"PRId64"  is <0 or >=%"PRId64") in im_get_pixel_as_double\n", y, im->ny);
        status->status = OUT_OF_BOUND;
        return val;
    }
    
    switch(im->image_type)
    {
        case ITYPE_FLOAT:
            val = GET_PIX(f, x, y);
            break;
            
        case ITYPE_DOUBLE:
            val = GET_PIX(d, x, y);
            break;
            
        case ITYPE_INT64:
            val = GET_PIX(i64, x, y);
            break;
            
        case ITYPE_INT32:
            val = GET_PIX(i32, x, y);
            break;
            
        case ITYPE_INT16:
            val = GET_PIX(i16, x, y);
            break;
            
        case ITYPE_INT8:
            val = GET_PIX(i8, x, y);
            break;
            
        case ITYPE_UINT64:
            val = GET_PIX(ui64, x, y);
            break;
            
        case ITYPE_UINT32:
            val = GET_PIX(ui32, x, y);
            break;
            
        case ITYPE_UINT16:
            val = GET_PIX(ui16, x, y);
            break;
            
        case ITYPE_UINT8:
            val = GET_PIX(ui8, x, y);
            break;
            
        default:
            fprintf(stderr, "Unknown/Invalid type passed to im_get_pixel_as_double (%d)\n", im->image_type);
            break;
    }
    
    return val;
    
}

