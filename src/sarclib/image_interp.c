/***************************************************************************
 * 
 *           Module :  image_interp.c
 *          Program :  sarclib
 *       Created by :  Mark Ashe on 14/08/2007
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions needed to interpolate images
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

#include "image_interp.h"

int generate_axis(double *axis, double axis_start, double axis_stop, long int np, double *output_resolution)
{	// simple axis generator
    //
	double res;
	if (np > 0)
	{	res	= (axis_stop - axis_start) / ( (double) (np-1) );
		long int i;
		for (i=0;i<np;i++) axis[i] = axis_start + (( (double)i ) * res);
		*output_resolution = res;
		return 1;
	}
	else return 0;
};

void simple_tpose(double *in, long int nx, long int ny, double *out)
{	// just swap the storage order 
    //
	long int i,j;
	for (i=0;i<nx;i++)
	{	for (j=0;j<ny;j++) out[j + i*ny] = in[i + nx*j];	}
};

long int index_below(double value, double *axis, long int np)
{	long int ind;
	if ((axis[np-1] - axis[0]) > 0)
	{	// data is forward sorted 
        //
        for (ind=1;ind<(np-1);ind++) if (axis[ind] >= value) return (ind-1);
	}
	else
	{	// data is reverse sorted 
        //
        for (ind=(np-2);ind>0;ind--) if (axis[ind] >= value) return (ind+1);
	}
	return ind;
};

long int index_above(double value, double *axis, long int np)
{	long int ind;
	if ((axis[np-1] - axis[0]) > 0)
	{	// data is forward sorted 
        //
        for (ind=(np-2);ind>0;ind--) if (axis[ind] <= value) return (ind+1);
	}
	else
	{	// data is reverse sorted 
        //
        for (ind=1;ind<(np-1);ind++) if (axis[ind] <= value) return (ind-1);
	}
	return ind;
};

int index_span(double minv, double maxv, double *axis, long int np, long int *below, long int *above)
{	// use previous functions until I get smarter... 
    //
    long int x0, x1;
	if (axis[np-1] > axis[0])
	{	x0	= index_below(minv, axis, np);
		x1	= index_above(maxv, axis, np);
	}
	else
	{	x1	= index_above(minv, axis, np);
		x0	= index_below(maxv, axis, np);
	}
	if (x0 <= x1)
	{	*below	= x0;
		*above	= x1;
	}
	else
	{	*below	= x1;
		*above	= x0;
	}
	return 1;
};

int index_bracketing(double value, double *axis_in, long int axis_length, long int *low_bound, long int *high_bound)
{	double dir;
	long int i;
	dir	= axis_in[1] - axis_in[0];
	if (dir > 0)
	{	// for forward sorted data 
        //
        if (value < axis_in[0])
		{	*low_bound	= 0;
			*high_bound	= 1;
		}
		else
		{	*low_bound	= axis_length-2;
			*high_bound	= axis_length-1;
			for (i=1; i<axis_length; i++)
			{	if ( (axis_in[i-1] < value)&&(axis_in[i] >= value) )
            {	*low_bound	= i-1;
                *high_bound	= i;
				return 1;
            }
			}
		}
	}
	else
	{	// for backward sorted data 
        //
        if (value > axis_in[0])
		{	*low_bound	= 0;
			*high_bound	= 1;
		}
		else
		{	*low_bound	= axis_length-2;
			*high_bound	= axis_length-1;
			for (i=1; i<axis_length; i++)
			{	if ( (axis_in[i-1] > value)&&(axis_in[i] <= value) )
            {	*low_bound	= i;
                *high_bound	= i-1;
                return 1;
            }
			}
		}
	}
	return 0;
};

int spline(double *in_x, double *in_y, long int np, double *out_y2)
{	// computes the spline stuff (deep magic) 
    //
    long int i;
	double u[np];
	double sig, p, qn, un;
	// making the assumption that this is a natural spline (U/K boundary y'') 
    //
    out_y2[0]	= 0;
	u[0]		= 0;
	qn		= 0;
	un		= 0;
	for (i=1;i<(np-1);i++)
	{	sig		= (in_x[i]-in_x[i-1])/(in_x[i+1]-in_x[i-1]);
		p		= sig*out_y2[i-1]+2.0;
		out_y2[i]	= (sig-1.0)/p;
		u[i]		= (in_y[i+1]-in_y[i])/(in_x[i+1]-in_x[i]);
		u[i]		= u[i] - (in_y[i]-in_y[i-1])/(in_x[i]-in_x[i-1]);
		u[i]		= (6.0*u[i]/(in_x[i+1]-in_x[i-1]) - sig*u[i-1])/p;
	}
	out_y2[np-1]	= (un-qn*u[(np-2)])/(qn*out_y2[np-1]+1.0);
	for (i=(np-2);i>=0;i--) out_y2[i] = out_y2[i]*out_y2[i+1] + u[i];
	return 1;
};

int reg_spline(double *in_x, double *in_y, long int np, double *out_y2)
{	// computes the spline stuff (deep magic) 
    //
    long int i;
	double u[np];
	double p, spc;
	
	// this routine is for REGULARLY GRIDDED DATA only 
    //
    spc	= (in_y[np-1] - in_y[0]) / ((double) np);
	
	// making the assumption that this is a natural spline (U/K boundary y'') 
    //
    out_y2[0]	= 0;
	u[0]		= 0;
    
	for (i=1;i<(np-1);i++)
	{	p		= 0.5*out_y2[i-1]+2.0;
		out_y2[i]	= -0.5/p;
		u[i]		= (in_y[i+1]-in_y[i])/spc;
		u[i]		= u[i] - (in_y[i]-in_y[i-1])/spc;
		u[i]		= (3.0*u[i]/spc - 0.5*u[i-1])/p;
	}
	out_y2[np-1]	= 0.0;
	for (i=(np-2);i>=0;i--) out_y2[i] = out_y2[i]*out_y2[i+1] + u[i];
	return 1;
};
int reg_splint(double *out_x, double *out_y,long int np, double x0, double x1, double y0, double y1, double dy0, double dy1)
{	// takes a set of arguments and splints accordingly 
    //
    long int i;
	double a, b, spc, h;
	
	spc	= x1 - x0;
	h	= spc*spc / 6.0;
	
	for (i=0;i<np;i++)
	{	a	= (x1-out_x[i])/spc;
		b	= (out_x[i]-x0)/spc;
		out_y[i] = a*y0 + b*y1 + h*((a*a*a-a)*dy0 + (b*b*b-b)*dy1);
	}
	return 0;
};

double splint(double out_x, double *in_x, double *in_y, long int np, double *in_y2)
{	// for a given out_x and spline construct in_y2 gives y 
    //
    long int klo,khi;
	double h,b,a, y_o;
    
	index_bracketing(out_x, in_x, np, &klo, &khi);
	
	h	= in_x[khi]-in_x[klo];
	if (h == 0.0) return -999;
	a	= (in_x[khi]-out_x)/h;
	b	= (out_x-in_x[klo])/h;
	y_o	= a*in_y[klo]+b*in_y[khi]+((a*a*a-a)*in_y2[klo]+(b*b*b-b)*in_y2[khi])*(h*h)/6.0;
	return y_o;
};

double blint_dbl(double x_out, double y_out, double x0, double x1, double y0, double y1, double z00, double z01, double z10, double z11)
{	// simplest possible billinear interpolation routine! 
	// for the inputs: you have x0, x1, y0, y1 that form a square so: 
	//	(x0, y1, z01)  (x1,y1, z11) 
	//	(x0, y0, z00)  (x1,y0, z10) 
    //
    double d1, d2, t, u, val;
	d1	= x1 - x0;
	d2	= y1 - y0;
	t	= (y_out - y0) / d2;
	u	= (x_out - x0) / d1;
	// this algorithm uses a grid square labelled so: 
	//	2 3	
	//	1 4	
    //
    val	= (1.0-t)*(1.0-u)*z00 + t*(1.0-u)*z01 + t*u*z11 + (1.0-t)*u*z10;
	return val;
};

int bilinterp_image(SPImage *inp_img, SPImage *ret_img, double *in_ties, double *out_ties, long int out_nx, long int out_ny, SPStatus *status)
{	// Demands rigorously formatted inputs:					
	//	inp_img and ret_img are VECTOR type with [y*nx + x] indexing	
	// 	in_ties and out_ties are [lon0, lat0, lon1, lat1]		
    //
	// the variables required for all interpolators				
    //
    long int x0,x1,y0,y1,x,y,ox0,ox1,oy0,oy1,ox,oy,indy0,indy1,in_nx, in_ny;
	double in_x[(*inp_img).nx], in_y[(*inp_img).ny], out_x[out_nx], out_y[out_ny], xres, yres;
	SPVector pxl;
	// Bilinear interpolation variables 
    //
    double d1,d2,z0,z1,z2,z3,u,t,val;
	
	// we use these a lot - make a local copy				
    //
    in_nx	= (*inp_img).nx;
	in_ny	= (*inp_img).ny;
	
	// set up the axis to interpolate from and to 				
    //
    generate_axis(in_x, in_ties[0], in_ties[2], in_nx, &xres);
	generate_axis(in_y, in_ties[1], in_ties[3], in_ny, &yres);
	
	generate_axis(out_x, out_ties[0], out_ties[2], out_nx, &xres);
	generate_axis(out_y, out_ties[1], out_ties[3], out_ny, &yres);
	
	// make sure the input image has the right resolution values 		
    //
    (*inp_img).xspc	= xres;
	(*inp_img).yspc	= yres;
	
	// the input array indices between which to work			
    //
    index_span(out_ties[0], out_ties[2], in_x, in_nx, &x0, &x1);
	index_span(out_ties[1], out_ties[3], in_y, in_ny, &y0, &y1);
	if (x0<1) x0++;
	if (y0<1) y0++;
	
	// BILINEAR INTERPOLATION 
    //
    for (y=y0;y<y1;y++)
	{	index_span(in_y[y-1], in_y[y], out_y, out_ny, &oy0, &oy1);
		d2	= in_y[y] - in_y[y-1];
		indy0	= (y-1)*in_nx;
		indy1	= indy0 + in_nx;
		for (x=x0;x<x1;x++)
		{	index_span(in_x[x-1], in_x[x], out_x, out_nx, &ox0, &ox1);
			d1	= in_x[x] - in_x[x-1];
			
			z0	= ((*inp_img).data.vect[indy0 + (x-1)]).z;
			z1	= ((*inp_img).data.vect[indy0 + x]).z;
			z2	= ((*inp_img).data.vect[indy1 + x]).z;
			z3	= ((*inp_img).data.vect[indy1 + (x-1)]).z;
            
			// only run over the output indices inside the span of the input ones 
            //
            for (oy=oy0;oy<=oy1;oy++)
			{	u	= (out_y[oy] - in_y[y-1]) / d2;
				for (ox=ox0;ox<=ox1;ox++)
				{	t	= (out_x[ox] - in_x[x-1]) / d1;
					val	= (1.0-t)*(1.0-u)*z0 + t*(1.0-u)*z1 + t*u*z2 + (1.0-t)*u*z3;
					VECT_CREATE(out_x[ox], out_y[oy], val, pxl);
					(*ret_img).data.vect[ox + oy*out_nx] = pxl;
				}
			}
		}
	}
	return 0;
};

int splinterp_image(SPImage *inp_img, SPImage *ret_img, double *in_ties, double *out_ties, long int out_nx, long int out_ny, SPStatus *status)
{	// Demands rigorously formatted inputs:					
	//	inp_img and ret_img are VECTOR type with [y*nx + x] indexing	
	// 	in_ties and out_ties are [lon0, lat0, lon1, lat1]		
    //
	// the variables required for the spline interpolator			
    //
    long int in_nx, in_ny, x, y, ox1, ox0;
	double in_x[(*inp_img).nx], in_y[(*inp_img).ny], out_x[out_nx], out_y[out_ny], xres, yres, val;
	SPVector pxl;
	int st;
	double x_spl[out_nx * ((*inp_img).ny)], dx2[(*inp_img).nx], dy2[(*inp_img).ny], tx_spl[out_nx * ((*inp_img).ny)], raw[((*inp_img).nx) * ((*inp_img).ny)];
	
	// we use these a lot - make a local copy				
    //
    in_nx	= (*inp_img).nx;
	in_ny	= (*inp_img).ny;
	
	// accept differnt types of data 
    //
    switch((*inp_img).image_type){
        case ITYPE_VECTOR :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = ((*inp_img).data.vect[y*in_nx + x]).z;
            }
            break;
		case ITYPE_UINT64 :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.ui64[y*in_nx + x]);
            }
            break;
		case ITYPE_UINT32 :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.ui32[y*in_nx + x]);
            }
            break;
		case ITYPE_UINT16 :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.ui16[y*in_nx + x]);
            }
            break;
		case ITYPE_UINT8  :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.ui8[y*in_nx + x]);
            }
            break;
		case ITYPE_INT64  :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.i64[y*in_nx + x]);
            }
            break;
		case ITYPE_INT32  :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.i32[y*in_nx + x]);
            }
            break;
		case ITYPE_INT16  :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.i16[y*in_nx + x]);
            }
            break;
		case ITYPE_INT8   :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.i8[y*in_nx + x]);
            }
            break;
		case ITYPE_FLOAT  :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (double)((*inp_img).data.f[y*in_nx + x]);
            }
            break;
		case ITYPE_DOUBLE :
            for (y=0;y<in_ny;y++){
                for (x=0;x<in_nx;x++)	raw[y*in_nx + x] = (*inp_img).data.d[y*in_nx + x];
            }
            break;
        default:
            printf("Sorry, type %d (%s) is not supported\n", inp_img->image_type, itype2string(inp_img->image_type) );
            status->status = INVALID_TYPE;
            return 1;
	}
	
	// set up the axis to interpolate from and to 				
    //
    generate_axis(in_x, in_ties[0], in_ties[2], in_nx, &xres);
	generate_axis(in_y, in_ties[1], in_ties[3], in_ny, &yres);
	generate_axis(out_x, out_ties[0], out_ties[2], out_nx, &xres);
	generate_axis(out_y, out_ties[1], out_ties[3], out_ny, &yres);
	(*ret_img).xspc	= xres;
	(*ret_img).yspc	= yres;
	
	for (y=0;y<in_ny;y++)
	{	st	= reg_spline(in_x, &(raw[y*in_nx]), in_nx, dx2);
		for (x=1;x<in_nx;x++)
		{	index_span(in_x[x-1], in_x[x], out_x, out_nx, &ox0, &ox1);
			reg_splint(&(out_x[ox0]), &(x_spl[y*out_nx + ox0]),(ox1-ox0), in_x[x-1], in_x[x], raw[y*in_nx+x-1], raw[y*in_nx+x], dx2[x-1], dx2[x]);
		}
	}
	
	simple_tpose(x_spl, out_nx, in_ny, tx_spl);
    
	// the Y-spline - column by column spline	
    //
    for (x=0;x<out_nx;x++)
	{	st	= spline(in_y, &(tx_spl[x*in_ny]), in_ny, dy2);
		for (y=0;y<out_ny;y++)
		{	val = splint(out_y[y], in_y, &(tx_spl[x*in_ny]), in_ny, dy2);
			VECT_CREATE(out_x[x], out_y[y], val, pxl);
			(*ret_img).data.vect[y*out_nx + x] = pxl;
		}
	}
	
	return 0;
};

int bilinterp_vecim(SPImage *inp_img, SPImage *ret_img, double *out_ties, long int out_nx, long int out_ny, SPStatus *status)
{	// only works on VECTOR type images					
    //

	// the variables required for all interpolators
    //
    long int x,y,ox0,ox1,oy0,oy1,ox,oy,indy0,indy1,in_nx, in_ny;
	double out_x[out_nx], out_y[out_ny], xres, yres, miny, maxy, minx, maxx;
	SPVector pxl;
	// Bilinear interpolation variables 
    //
    double d1,d2,z0,z1,z2,z3,t,u,val;
	double x0,x1,y0,y1;
	
	// we use these a lot - make a local copy				
    //
    in_nx	= (*inp_img).nx;
	in_ny	= (*inp_img).ny;
	
	// the input array indices between which to work cannot easily be 	
	// defined due to the potential for irregular gridding			
    //
    generate_axis(out_x, out_ties[0], out_ties[2], out_nx, &xres);
	generate_axis(out_y, out_ties[1], out_ties[3], out_ny, &yres);
	
	// the image to be filled						
    //
    im_create(ret_img, ITYPE_VECTOR, out_nx, out_ny, xres, yres, status);
	
	// get a representative 'xres' and 'yres'				
    //
    x	= in_nx / 2;
	y	= in_ny / 2;
	y0	= ((*inp_img).data.vect[y*in_nx + x]).y - ((*inp_img).data.vect[(y-1)*in_nx + x]).y;
	x0	= ((*inp_img).data.vect[y*in_nx + x]).x - ((*inp_img).data.vect[y*in_nx + x-1]).x;
	y0	= fabs(y0);
	x0	= fabs(x0);
	
	if(out_ties[1] < out_ties[3])
	{	miny	= out_ties[1];
		maxy	= out_ties[3];
	}
	else
	{	miny	= out_ties[3];
		maxy	= out_ties[1];
	}
	if(out_ties[0] < out_ties[2])
	{	minx	= out_ties[0];
		maxx	= out_ties[2];
	}
	else
	{	minx	= out_ties[2];
		maxx	= out_ties[0];
	}
	miny	= miny - 2.0*y0;
	maxy	= maxy + 2.0*y0;
	minx	= minx - 2.0*x0;
	maxx	= maxx + 2.0*x0;
	
	// BILINEAR INTERPOLATION 
    //
    for (y=1;y<in_ny;y++)
	{	indy0	= (y-1)*in_nx;
		indy1	= indy0 + in_nx;
		y0	= ((*inp_img).data.vect[indy0]).y;
		y1	= ((*inp_img).data.vect[indy1]).y;
		
		// check BOTH points are within bounds			
		// this check needs work - last point may not quite be in bounds 
        //
        if (( y0 <= maxy )&&( y0 >= miny )&&( y1 <= maxy )&&( y1 >= miny ))
		{	d2	= y1-y0;
			index_span(y0, y1, out_y, out_ny, &oy0, &oy1);
            
			for (x=1;x<in_nx;x++)
			{	x0	= ((*inp_img).data.vect[indy0+x-1]).x;
				x1	= ((*inp_img).data.vect[indy1+x]).x;
                
				// check BOTH points are within bounds			
                //
                if ((x0 <= maxx )&&( x0 >= minx )&&( x1 <= maxx )&&( x1 >= minx))
				{	d1	= x1-x0;
					z0	= ((*inp_img).data.vect[indy0 + (x-1)]).z;
					z1	= ((*inp_img).data.vect[indy0 + x]).z;
					z2	= ((*inp_img).data.vect[indy1 + x]).z;
					z3	= ((*inp_img).data.vect[indy1 + (x-1)]).z;
					
					index_span(x0, x1, out_x, out_nx, &ox0, &ox1);
                    
					// only run over the output indices inside the span of the input ones 
                    //
                    for (oy=oy0;oy<=oy1;oy++)
					{	u	= (out_y[oy] - y0) / d2;
						for (ox=ox0;ox<=ox1;ox++)
						{	t	= (out_x[ox] - x0) / d1;
							val	= (1.0-t)*(1.0-u)*z0 + t*(1.0-u)*z1 + t*u*z2 + (1.0-t)*u*z3;
							VECT_CREATE(out_x[ox], out_y[oy], val, pxl);
							(*ret_img).data.vect[ox + oy*out_nx] = pxl;
						}
					}
				}
			}
		}
	}
	return 0;
};

int splinterp_vecim(SPImage *inp_img, SPImage *ret_img, double *out_ties, long int out_nx, long int out_ny, SPStatus *status)
{	// only works on VECTOR type images					
    //

	// the variables required for spline interpolators
    //
    long int x,y,in_nx,in_ny,y_off,y_off2,ox0,ox,ox1;
	double out_x[out_nx], out_y[out_ny], xres, yres, miny, maxy, minx, maxx;
	SPVector pxl;
	double val,sig,p,h,a,b,x_1,x0,x1,y_1,y0,y1,d0,d_1,y_3d;
	double x_u[(*inp_img).nx], dy2[(*inp_img).nx], dx2[(*inp_img).ny], y_u[(*inp_img).ny];
	SPImage spl_x;
	
	// we use these a lot - make a local copy				
    //
    in_nx	= (*inp_img).nx;
	in_ny	= (*inp_img).ny;
	
	// the input array indices between which to work cannot easily be 	
	// defined due to the potential for irregular gridding			
    //
    generate_axis(out_x, out_ties[0], out_ties[2], out_nx, &xres);
	generate_axis(out_y, out_ties[1], out_ties[3], out_ny, &yres);
    
	im_create(&spl_x, ITYPE_VECTOR, out_nx, in_ny, 1.0, 1.0, status);
	im_create(ret_img, ITYPE_VECTOR, out_nx, out_ny, xres, yres, status);
	
	// get a representative 'xres' and 'yres'				
    //
    x	= in_nx / 2;
	y	= in_ny / 2;
	y0	= ((*inp_img).data.vect[y*in_nx + x]).y - ((*inp_img).data.vect[(y-1)*in_nx + x]).y;
	x0	= ((*inp_img).data.vect[y*in_nx + x]).x - ((*inp_img).data.vect[y*in_nx + x-1]).x;
	y0	= fabs(y0);
	x0	= fabs(x0);
	
	if(out_ties[1] < out_ties[3])
	{	miny	= out_ties[1];
		maxy	= out_ties[3];
	}
	else
	{	miny	= out_ties[3];
		maxy	= out_ties[1];
	}
	if(out_ties[0] < out_ties[2])
	{	minx	= out_ties[0];
		maxx	= out_ties[2];
	}
	else
	{	minx	= out_ties[2];
		maxx	= out_ties[0];
	}
	miny	= miny - 2.0*y0;
	maxy	= maxy + 2.0*y0;
	minx	= minx - 2.0*x0;
	maxx	= maxx + 2.0*x0;
	
	// SPLINE INTERPOLATION 
    //
    for (y=0;y<in_ny;y++)
	{	// for every line....								
        //
        y_off	= y * in_nx;
		y_off2	= y * out_nx;
		// get a dy2/dx2 (assumption: natural spline, dy2dx2 = 0 at boundaries		
        //
        dy2[0]	= 0.00;
		x_u[0]	= 0.00;
		dy2[in_nx-1] = 0.00;
		for (x=1;x<(in_nx-1);x++)
		{	// the x-axis of a 2D spline (= data x axis)				
            //
            x_1	= ((*inp_img).data.vect[y_off + x-1]).x;
			x0	= ((*inp_img).data.vect[y_off + x]).x;
			x1	= ((*inp_img).data.vect[y_off + x+1]).x;
			// the y-axis of a 2D spline (= data values / z-axis)			
            //
            y_1	= ((*inp_img).data.vect[y_off + x-1]).z;
			y0	= ((*inp_img).data.vect[y_off + x]).z;
			y1	= ((*inp_img).data.vect[y_off + x+1]).z;
			sig	= (x0 - x_1) / (x1 - x_1);
			p	= sig * dy2[x-1] + 2.0;
			dy2[x]	= (sig - 1.0) / p;
			x_u[x]	= ((y1 - y0) / (x1 - x0)) - ((y0 - y_1) / (x0 - x_1));
			x_u[x]	= (6.0*x_u[x] / (x1-  x_1)) - (sig * x_u[x-1] / p);
		}
		// forward loop to get the 'true' dy2dx2					
        //
        for (x = (in_nx-2); x>=0; x--) dy2[x] = (dy2[x] * dy2[x+1]) + x_u[x];
		// now use this dy2 to get a splinted value for the axis			
		// by looping over every input axis point...					
        //
        for (x=1;x<in_nx;x++)
		{	x_1	= ((*inp_img).data.vect[y_off + x-1]).x;
			x0	= ((*inp_img).data.vect[y_off + x]).x;
            
			// to activate these should both be inside interpolation bounds 	
            //
            if ( (x_1>=minx)&&(x_1<=maxx)&&(x0>=minx)&&(x0<=maxx) )
			{	// finish allocating temporary variables			
                //
                y_1	= ((*inp_img).data.vect[y_off + x-1]).z;
				y0	= ((*inp_img).data.vect[y_off + x]).z;
				d_1	= dy2[x-1];
				d0	= dy2[x];
				// set up a points true y value (using the mean of the two ys)	
                //
                y_3d	=  ((*inp_img).data.vect[y_off + x-1]).y;
				y_3d	=  y_3d + ((*inp_img).data.vect[y_off + x]).y;
				y_3d	= y_3d / 2.0;
				// find members of output axis that are inside the x values 	
                //
                index_span(x_1, x0, out_x, out_nx, &ox0, &ox1);
				h	= x0 - x_1;
				// loop over these and splint					
                //
                for (ox=ox0;ox<=ox1;ox++)
				{	a	= (x0 - out_x[ox]) / h;
					b	= (out_x[ox] - x_1) / h;
					val	= (a * y_1) + (b * y0) + ((((a*a*a-a)*d_1) + ((b*b*b-b)*d0))*h*h/6.0);
					VECT_CREATE(out_x[ox], y_3d, val, pxl);
					spl_x.data.vect[y_off2 + ox] = pxl;
				}
			}
		}
	}
	
	// now we have to spline along the OTHER axis						
    //
    im_transp (&spl_x, status);
	
	// there are now out_nx rows of in_ny points....					
    //
    for (y=0;y<out_nx;y++)
	{	// for every line....								
        //
        y_off	= y * in_ny;
		y_off2	= y * out_ny;
		// get a dx2/dy2 (assumption: natural spline, dx2dy2 = 0 at boundaries		
        //
        dx2[0]	= 0.00;
		y_u[0]	= 0.00;
		dx2[in_ny-1] = 0.00;
		for (x=1;x<(in_ny-1);x++)
		{	// the x-axis of a 2D spline (= interpolated values / y axis)		
            //
            x_1	= (spl_x.data.vect[y_off + x-1]).y;
			x0	= (spl_x.data.vect[y_off + x]).y;
			x1	= (spl_x.data.vect[y_off + x+1]).y;
			// the y-axis of a 2D spline (= interpolated values / z-axis)		
            //
            y_1	= (spl_x.data.vect[y_off + x-1]).z;
			y0	= (spl_x.data.vect[y_off + x]).z;
			y1	= (spl_x.data.vect[y_off + x+1]).z;
			sig	= (x0 - x_1) / (x1 - x_1);
			p	= sig * dx2[x-1] + 2.0;
			dx2[x]	= (sig - 1.0) / p;
			y_u[x]	= ((y1 - y0) / (x1 - x0)) - ((y0 - y_1) / (x0 - x_1));
			y_u[x]	= (6.0*y_u[x] / (x1-  x_1)) - (sig * y_u[x-1] / p);
		}
		// forward loop to get the 'true' dx2dy2					
        //
        for (x = (in_ny-2); x>=0; x--) dx2[x] = (dx2[x] * dx2[x+1]) + y_u[x];
		// now use this dy2 to get a splinted value for the axis			
		// by looping over every input axis point...					
        //
        for (x=1;x<in_ny;x++)
		{	x_1	= (spl_x.data.vect[y_off + x-1]).y;
			x0	= (spl_x.data.vect[y_off + x]).y;
            
			// to activate these should both be inside interpolation bounds 	
            //
            if ( (x_1>=miny)&&(x_1<=maxy)&&(x0>=miny)&&(x0<=maxy) )
			{	// finish allocating temporary variables			
                //
                y_1	= (spl_x.data.vect[y_off + x-1]).z;
				y0	= (spl_x.data.vect[y_off + x]).z;
				d_1	= dx2[x-1];
				d0	= dx2[x];
				
				// set up a points true y value (using the mean of the two xs)	
                //
                y_3d	=  (spl_x.data.vect[y_off + x-1]).x;
				y_3d	=  y_3d + (spl_x.data.vect[y_off + x-1]).x;
				y_3d	= y_3d / 2.0;
				// find members of output axis that are inside the x values 	
                //
                index_span(x_1, x0, out_y, out_ny, &ox0, &ox1);
				h	= x0 - x_1;
				// loop over these and splint					
                //
                for (ox=ox0;ox<=ox1;ox++)
				{	a	= (x0 - out_y[ox]) / h;
					b	= (out_y[ox] - x_1) / h;
					val	= (a * y_1) + (b * y0) + ((((a*a*a-a)*d_1) + ((b*b*b-b)*d0))*h*h/6.0);
					VECT_CREATE(out_y[ox], y_3d, val, pxl);
					(*ret_img).data.vect[ox*out_ny + y] = pxl;
				}
			}
		}
	}
    
	im_destroy(&spl_x,status);
	return 0;
};

