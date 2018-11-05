/***************************************************************************
 * 
 *           Module :  interpolate.c
 *          Program :  sarclib
 *       Created by :  Alan Blake on 23/08/2005
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *       Sinc interpolation routines.
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

#include "sarclib.h"

static void ham1dx(double * data, int nx);
static int ilib_nint(double x);



/***************************************************************************
 * 
 *           Module :  interpolate.c
 *          Program :  sarclib
 *       Created by :  Alan Blake on 23/08/2005
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *       Sinc interpolation routines.
 *
 * Description:   Uses a 8 point interpolation
 *                function.
 *
 *****************************************************/
int
interpolate8pt(SPImage * data, SPImage * new_image)
{
    int x;
    int kernel_nx;
    
    double k_start, k_end, k_mid;
    double * kernel;
    double xpf, k_xpf;
    double xps, xpe;
    int xp;
    int oversample_factor = 100;
    double k_osf = 1.0 / (double) oversample_factor;
    int npoint = 32;
    
    kernel_nx = oversample_factor * npoint + 1;
    kernel = generate_sinc_kernel(oversample_factor, npoint);
    
    k_start = 1.0 - (kernel_nx / 2) * k_osf;
    k_end = 1.0 + ((kernel_nx - 1) - (kernel_nx / 2)) * k_osf;
    k_mid = (0.5 + (kernel_nx / 2)) * k_osf;
    
    for(x = 0; x < new_image->nx; x++)
    {
        xpf = (double) x * (double) data->nx / (double) new_image->nx;
        
        xps = (int)(xpf + k_start);
        if (xps < 0) xps = 0;
        
        xpe = (int)(xpf + k_end);
        if (xpe > data->nx) xpe = data->nx;
        
        new_image->data.cmpl_f[x].r = 0.0;
        new_image->data.cmpl_f[x].i = 0.0;
        k_xpf = (xps - xpf + k_mid) / k_osf;
        
        for(xp = xps; xp < xpe; xp++)
        {
            new_image->data.cmpl_f[x].r += data->data.cmpl_f[xp].r * kernel[(int)k_xpf];
            new_image->data.cmpl_f[x].i += data->data.cmpl_f[xp].i * kernel[(int)k_xpf];
            k_xpf +=  1.0 / k_osf;
        }
        
    }
    free(kernel);
    
    return(0);
}

/***********************************************************
 *
 * This routine is used to interpolate a supplied image
 * using a kernel function supplied as an oversampled image
 * to give an output image where the values correspond with
 * positions supplied in a complex position image. Interpolation
 * can be restricted to nearest neighbour in either or both
 * directions. If kernel interpolation in x and y is required then
 * the oversampled kernel can be a one or full 2d image.
 *
 */
typedef enum {
	NX_NY,   /* nearest neighbour in x , nearest neighbour in y */
	NX_KY,   /* nearest neighbour in x , kernel in y            */
	KX_NY,   /* kernel in x            , nearest neighbour in y */
	KX_KY    /* kernel in x            , kernel in y            */
} INTERP_MODE;

#define KX_KY_MACRO(clear_sum,active_line)                                  \
{                                                                           \
    long int xp, xps, xpe;                                                  \
    float k_xpf, k_xpf_tmp;                                                 \
                                                                            \
    long int yp, yps, ype;                                                  \
    float k_ypf;                                                            \
                                                                            \
    double ky;                                                              \
                                                                            \
    for (j=0; j < nj ; j++) {                                               \
        for(i=0; i < ni ; i++) {                                            \
                                                                            \
            xpf=in_rsx*m_data[i].r;                                         \
            xps=(int)(xpf + k_start); if(xps < 0)       xps = 0;            \
            xpe=(int)(xpf + k_end);   if(xpe > in_xlim) xpe = in_xlim;      \
                                                                            \
            ypf=in_rsy*m_data[i].i;                                         \
            yps=(int)(ypf + k_start); if(yps < 0)       yps = 0;            \
            ype=(int)(ypf + k_end );  if(ype > in_ylim) ype = in_ylim;      \
                                                                            \
            k_xpf=(xps-xpf+k_mid)*k_rs;                                     \
            k_ypf=(yps-ypf+k_mid)*k_rs;                                     \
                                                                            \
            in_data=data->data.cmpl_f+yps*in_stride;                        \
                                                                            \
            {clear_sum}                                                     \
                                                                            \
            for(yp=yps; yp < ype; yp++) {                                   \
                k_xpf_tmp=k_xpf;                                            \
                ky=k_data[(int)k_ypf];                                      \
                for(xp=xps; xp < xpe; xp++) {                               \
                    {active_line}                                           \
                    k_xpf_tmp+=k_rs;                                        \
                }                                                           \
                k_ypf+=k_rs;                                                \
                in_data+=in_stride;                                         \
            }                                                               \
            out_data[i]=sum;                                                \
        }                                                                   \
        out_data+=out_stride;                                               \
        m_data+=m_stride;                                                   \
    }                                                                       \
}

#define KX_NY_MACRO(clear_sum,clear_result,active_line)                     \
{                                                                           \
    long int xp, xps, xpe;                                                  \
    float k_xpf;                                                            \
                                                                            \
    for (j=0; j < nj ; j++) {                                               \
        for(i=0; i < ni ; i++) {                                            \
            yp=ilib_nint(in_rsy*m_data[i].i);                               \
                                                                            \
            if(yp < 0 || yp >= in_ylim ) {clear_result} else {              \
                                                                            \
                xpf=in_rsx*m_data[i].r;                                     \
                xps=(int)(xpf + k_start); if(xps < 0)       xps = 0;        \
                xpe=(int)(xpf + k_end);   if(xpe > in_xlim) xpe = in_xlim;  \
                                                                            \
                {clear_sum}                                                 \
                                                                            \
                k_xpf=(xps-xpf+k_mid)*k_rs;                                 \
                in_data=data->data.cmpl_f+yp*in_stride;                     \
                                                                            \
                for(xp=xps; xp < xpe; xp++) {                               \
                    {active_line}                                           \
                    k_xpf+=k_rs;                                            \
                }                                                           \
                out_data[i]=sum;                                            \
            }                                                               \
        }                                                                   \
        out_data+=out_stride;                                               \
        m_data+=m_stride;                                                   \
    }                                                                       \
}


#define NX_KY_MACRO(clear_sum,clear_result,active_line)                     \
{                                                                           \
    long int yp, yps, ype;                                                  \
    float k_ypf;                                                            \
                                                                            \
    for (j=0; j < nj ; j++) {                                               \
        for(i=0; i < ni ; i++) {                                            \
            xp=ilib_nint(in_rsx*m_data[i].r);                               \
                                                                            \
            if(xp < 0 || xp >= in_xlim ) {clear_result} else {              \
                                                                            \
                ypf=in_rsy*m_data[i].i;                                     \
                yps=(int)(ypf + k_start); if(yps < 0)       yps = 0;        \
                ype=(int)(ypf + k_end);   if(ype > in_ylim) ype = in_ylim;  \
                                                                            \
                {clear_sum}                                                 \
                                                                            \
                k_ypf=(yps-ypf+k_mid)*k_rs;                                 \
                in_data=data->data.cmpl_f+yps*in_stride+xp;                 \
                                                                            \
                for(yp=yps; yp < ype; yp++) {                               \
                    {active_line}                                           \
                    k_ypf+=k_rs;                                            \
                    in_data+=in_stride;                                     \
                }                                                           \
                out_data[i]=sum;                                            \
            }                                                               \
        }                                                                   \
        out_data+=out_stride;                                               \
        m_data+=m_stride;                                                   \
    }                                                                       \
}

#define NX_NY_MACRO(clear_result,active_line)                               \
{                                                                           \
    in_data = data->data.cmpl_f;                                            \
    for (j=0; j < nj ; j++) {                                               \
        for(i=0; i < ni ; i++) {                                            \
            xp=ilib_nint(in_rsx*m_data[i].r);                               \
            yp=ilib_nint(in_rsy*m_data[i].i);                               \
            if(xp < 0 || xp >= in_xlim || yp < 0 || yp >= in_ylim ) {       \
                {clear_result}                                              \
            } else {                                                        \
                {active_line}                                               \
            }                                                               \
        }                                                                   \
        out_data+=out_stride;                                               \
        m_data+=m_stride;                                                   \
    }                                                                       \
}

int
image_interpolate(SPImage * data, SPImage *map_data, double * k_data, int kernel_nx, int kernel_oversample,
                  int use_kernel_x, int use_kernel_y, SPImage *o_data)
{
	float in_rsx = 1.0, in_rsy = 1.0;                   /* input reciprocal sample spacings, used to calculate input pixel positions in m*/
	long int ni, nj;                                    /* mapping image and hence output image size */
	float k_rs = kernel_oversample;                     /* reciprocal of kernel sample spacing */
	float k_start, k_end;                               /* kernel start and end values based on k_sx and number of points in kernel data */
	float k_mid;                                        /* mid point of kernel with offsets so truncation to float gives correct index into kernel array */
	long int in_stride;                                 /* input array stride */
	long int out_stride;                                /* output array stride */
	long int m_stride;                                  /* mapping array stride */
	int i, j;                                           /* loop over mapping/ output image */
	float xpf, ypf;                                     /* position to interpolate to in input image */
	int xp, yp;                                         /* position to interpolate to in input image, normalised to nearest pixel coordinates */
	INTERP_MODE imode;                                  /* type of interpolation to perform */
	long int in_xlim = data->nx, in_ylim = data->ny;    /* limiting maximum pixel values for input image */
	SPCmplx sum;
	double k_osf = 1.0 / (double) kernel_oversample;
	SPCmplx * in_data = data->data.cmpl_f;
	SPCmplx * out_data = o_data->data.cmpl_f;
	SPCmplx * m_data = map_data->data.cmpl_f;
	
	k_start = 1.0 - (kernel_nx / 2) * k_osf;
	k_end = 1.0 + ((kernel_nx - 1) - (kernel_nx / 2)) * k_osf;
	k_mid = (0.5 + (kernel_nx / 2)) * k_osf;
	
	ni=map_data->nx;
	nj=map_data->ny;
	
	in_stride = data->nx;
	out_stride = map_data->nx;
	m_stride = map_data->nx;
	
	imode = use_kernel_x ? (use_kernel_y ? KX_KY : KX_NY) : (use_kernel_y ? NX_KY : NX_NY);
	
	switch (imode) {
        case NX_NY:
            NX_NY_MACRO(out_data[i].r=0.0f;out_data[i].i=0.0f;,
                        out_data[i]=*(in_data+yp*in_stride+xp);) break;
        case NX_KY:
            NX_KY_MACRO(sum.r=0.0f;sum.i=0.0f;,
                        out_data[i].r=0.0f;out_data[i].i=0.0f;,
                        sum.r+=in_data->r * k_data[(int)k_ypf];sum.i+=in_data->i*k_data[(int)k_ypf];) break;
        case KX_NY:
            KX_NY_MACRO(sum.r=0.0f;sum.i=0.0f;,
                        out_data[i].r=0.0f;out_data[i].i=0.0f;,
                        sum.r+=in_data[xp].r*k_data[(int)k_xpf];sum.i+=in_data[xp].i*k_data[(int)k_xpf];) break;
        case KX_KY:
            KX_KY_MACRO(sum.r=0.0f; sum.i=0.0f;,
                        sum.r+=in_data[xp].r*ky*k_data[(int)k_xpf_tmp];sum.i+=in_data[xp].i*ky*k_data[(int)k_xpf_tmp];) break;
        default:
            fprintf(stderr, "Invalid mode in image_interpolate\n");
            exit(89);
	}
    
	return (0);
    
}



/*****************************************************
 *
 * Function name: interpolate_rotate
 * Input:         a data array
 *                The current size of the data
 *                The requested rotation (radians)
 * Output:        None.
 * Return:        The interpolated array.
 *
 * Description:   Uses a 8 point interpolation
 *                function to perform the rotation.
 *
 *****************************************************/
SPCmplx *
interpolate_rotate(SPCmplx * data, int nx, int ny, double theta)
{
    int x, y;
    int kernel_nx;
    
    SPCmplx * new_image;
    double k_start, k_end, k_mid;
    double * kernel;
    double xpf, k_xpf, k_xpf_tmp;
    double xps, xpe;
    int xp;
    double yps, ype;
    int yp;
    double ky, k_ypf, ypf;
    
    int oversample_factor = 100;
    double k_osf = 1.0 / (double) oversample_factor;
    double cos_theta, sin_theta;
    
    
    new_image = calloc(nx * ny, sizeof(SPCmplx));
    if (!new_image)
    {
        fprintf(stderr, "Failed to calloc %s:%d\n", __FILE__, __LINE__);
        exit(801);
    }
    kernel_nx = oversample_factor * 8 + 1;
    kernel = generate_sinc_kernel(oversample_factor, 8);
    cos_theta = cos(M_PI/2.0 - theta);
    sin_theta = sin(M_PI/2.0 - theta);
    
    k_start = 1.0 - (kernel_nx / 2) * k_osf;
    k_end = 1.0 + ((kernel_nx - 1) - (kernel_nx / 2)) * k_osf;
    k_mid = (0.5 + (kernel_nx / 2)) * k_osf;
    
    for(y = 0; y < ny; y++)
    {
        for(x = 0; x < nx; x++)
        {
            xpf = ((double) (x - nx/2) * cos_theta + (y - ny/2) * sin_theta) + nx/2;
            ypf = ((double) -(y - ny/2) * sin_theta + (x - nx/2) * cos_theta) + ny/2;
            
            xps = (int)(xpf + k_start);
            if (xps < 0) xps = 0;
            
            xpe = (int)(xpf + k_end);
            if (xpe > nx) xpe = nx;
            
            yps = (int)(ypf + k_start);
            if (yps < 0) yps = 0;
            
            ype = (int)(ypf + k_end);
            if (ype > ny) ype = ny;
            
            new_image[x + y * nx].r = 0.0;
            new_image[x + y * nx].i = 0.0;
            
            
            k_xpf = (xps - xpf + k_mid) / k_osf;
            k_ypf = (yps - ypf + k_mid) / k_osf;
            
            for (yp = yps; yp < ype; yp++)
            {
                ky = kernel[(int)k_ypf];
                k_xpf_tmp =  k_xpf;
                
                for(xp = xps; xp < xpe; xp++)
                {
                    new_image[x + y * nx].r += data[(int)(xp + yps * nx)].r * ky * kernel[(int)k_xpf_tmp];
                    new_image[x + y * nx].i += data[(int)(xp + yps * nx)].i * ky * kernel[(int)k_xpf_tmp];
                    k_xpf_tmp +=  1.0 / k_osf;
                }
                k_ypf +=  1.0 / k_osf;
            }
        }
    }
    
    return(new_image);
    
}

/*****************************************************
 *
 * Function name: generate_sinc_kernel
 * Input:         the oversample factor
 *                the number of points
 * Output:        None.
 * Return:        The sinc kernel.
 *
 * Description:   Generates a sinc kernel to use in
 *                the interpolation stuff.
 *
 *****************************************************/
double *
generate_sinc_kernel(int oversample_factor, int num_of_points)
{
    int n = oversample_factor * num_of_points + 1;
    int x;
    double val;
    double * data;
    double k = M_PI / (double) oversample_factor;
    
    if (oversample_factor == 0 || num_of_points == 0)
    {
        fprintf(stderr, "Either oversample_factor (= %d) or number of points (= %d) is zero in generate_sinc_kernel\n",
                oversample_factor, num_of_points);
        exit(803);
    }
    
    data = calloc(n, sizeof(double));
    if (!data)
    {
        fprintf(stderr, "Failed to calloc %s:%d\n", __FILE__, __LINE__);
        exit(802);
    }
    
    for(x = 0; x < n; x++)
    {
        val = k * (double)(x - n/2);
        data[x] = (fabs(val) < 1.0e-4 ) ? 1.0 : sin(val)/val;
    }
    
    ham1dx(data, n);
    
    return(data);
}

/*****************************************************
 *
 * Function name: ham1dx
 * Input:         The data aray
 *                the number of points
 * Output:        the weighted data
 * Return:        None
 *
 * Description:   Takes an input data array and
 *                multiplied by a weighting function.
 *
 *****************************************************/
static void
ham1dx(double * data, int nx)
{
    int x;
    
    double ped = 0.08;
    double a = 1.0 - ped;
    double val;
    
    for(x = 0; x < nx; x++)
    {
        val = M_PI * (-0.5 + (double) x / (double) (nx - 1));
        data[x] *= ped + a * cos(val) * cos(val);
    }
}

static int
ilib_nint(double x)
{
	return ( x>0.0 ? x+0.5 : x-0.5);
}
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
