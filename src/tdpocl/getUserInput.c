/***************************************************************************
 * 
 *           Module :  getUserInput.c
 *          Program :  tdpocl
 *       Created by :  Darren Muff on 14/03/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *   Parse the user input and return it back in a sensible way to the main tdp routine
 *   Heres the logic for parsing the input:
 *      USEGPU ?
 *      Use Topotrim File?
 *      If yes:
 *          Get TTM File
 *          Process Mission or Reference Image ?
 *          get cphd fname, Surface_fname, Nx, Ny, xstart, ystart
 *          How many pulses to process?
 *          Start pulse index ?
 *          Output Image Name?
 *      else:
 *          Get CPHD File
 *          Use a surface file (or default focal plane for image)
 *          If yes:
 *              get surface filename
 *              check it
 *              Get Nx  (default is surface NX)
 *              Get Ny  (default is surface NY)
 *              Get xStart (default (surfaceX - NX)/2 )
 *              Get yStart (default (surfaceY - NY)/2 )
 *              How many pulses to process?
 *              Start pulse index ?
 *              Output Image Name?
 *          else:
 *              Get SRP (default from CPHD file)
 *              Get Nx
 *              Get Ny
 *              Get xSpacing (m)
 *              Get ySPacing (m)
 *              How many pulses to process?
 *              Start pulse index ?
 *              Name of a file to save the generated surface to
 *              Output Image Name? (here rather than end because we want all input b4 we build surface file (slow)
 *              Generate surface and save it out
 *              store surface name Nx, Ny, xstart=0, ystart=0
 *          end if
 *      end if
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

#include <stdio.h>
#include "tdp.h"

void getPulses(CPHDHeader *hdr, int *nPulses, int *startPulse){
    *nPulses = hdr->num_azi ;
    char prompt[255];
    sprintf(prompt,"Number of pulses to process (1-%d)",hdr->num_azi);
    *nPulses = input_int(prompt, "NumPulse",
                         "Number of pulses in synthetic aperture to process when forming image", *nPulses);
    if(*nPulses > hdr->num_azi) {
        printf("*** Warning : Too many pulses requested (%d). File only has %d pulses. setting nPulses to maximum\n",*nPulses,hdr->num_azi);
        *nPulses = hdr->num_azi ;
    }
    
    *startPulse = (hdr->num_azi- *nPulses ) / 2 ;
    *startPulse = input_int("Start pulse", "StartPulse",
                            "First pulse in synthetic aperture to process", (hdr->num_azi - *nPulses)/2);
    if(*startPulse < 0){
        printf("*** Warning. Start pulse is negative (%d). Setting to zero\n",*startPulse);
        *startPulse = 0;
    }
    if (*startPulse + *nPulses > hdr->num_azi){
        printf("*** Error: Number of pulses (%d) and start pulse (%d) is greater than the available pulses in the dataset (%d)\n",*nPulses, *startPulse, hdr->num_azi);
        exit(1);
    }

    return ;
}

void getUserInput(int *USEGPU,              // Whether to use the GPU if its available or run on CPU
                  char **cphdFilename,      // Name of CPHD file to process
                  char **surfaceFilename,   // Name of surface file - defining each pixel in ECEF coords
                  int *Nx,                  // Size of output image in X
                  int *Ny,                  // Size of output image in Y
                  int *xstart,              // First x sample in surface file to use
                  int *ystart,              // First y sample in surface file to use
                  int *nPulses,             // Number of pulses in CPHD file to process
                  int *startPulse,          // First pulse to proces
                  double *EIRP,             // If non zero then use this value of Pt*Gt to output image as RCS (m^2)
                  char **outputFilename,    // Name of file to store the output image
                  SPStatus *status )        // mtlib status flag
{
    int useSurface ;
    SPImage im;
    CPHDHeader hdr;
    SPVector sc_pos;
    SPVector sc_ll;
    double surface_sx, surface_sy; 
    double slant_rng_res ;
        
    // *USEGPU   = input_yesno("Use GPU?","USEGPU",
    //                        "USE GPU if its available. If it isnt then use the CPU.",
    //                        TRUE);
    
    *cphdFilename = input_string("CPHD Filename", "InputFName",
                                 "Name of a CPHD file to process into an image",
                                 "/local_storage/DGM/first.cph");
    // Check the integrity of the CPHD file
    //
    im_init(&im, status);
    readCPHDHeader(*cphdFilename, &hdr, status);
    load_cphd(&im, &hdr, 0, 1, status);
    im_destroy(&im, status);
    
    getPulses(&hdr, nPulses, startPulse);

    useSurface = input_yesno("Use a surface file (or default focal plane for image)", "UseSurfaceFile",
                             "Use a surface file already generated. If no then a flat focal plane will be created", 0);
    
    if(!useSurface){
        
        // Find out the scene centre (stabalisation point), convert to a lat/long and then use that
        // as the default to ask the user where they want the scene centre to be
        //
        ecef2latlon(&hdr.pulses[*startPulse + *nPulses/2].srp, &sc_pos, status);
        sc_ll = input_vect("Scene centre lat/long/alt", "Scene centre",
                           "Enter the latitude, longitude and altitude as 3 decimal numbers seperated by spaces", sc_pos);
        
        // To reduce errors caused by latlon2ecef translation, then only convert if required
        //
        if( !(sc_ll.x == sc_pos.x && sc_ll.y == sc_pos.y && sc_ll.z == sc_pos.z) ){
            latlon2ecef(&sc_ll, &sc_pos, status);
        }else{
            sc_pos = hdr.pulses[*startPulse + *nPulses/2].srp ;
        }
        
        // Ask the user for the other parameters of the final image
        //
        slant_rng_res = floor((SIPC_c / (2.0 * hdr.pulse_length * hdr.chirp_gamma)) * 85.0) / 100.0;
        
        surface_sx = input_dbl("Image spacing width (m)", "IMsx",
                               "Cross range image spacing in metres", slant_rng_res);
        surface_sy = input_dbl("Image spacing height (m)", "IMsy",
                               "Range image spacing in metres", slant_rng_res);
        
        *Nx = input_int("Image width (pixels)", "IMnx",
                        "The cross range dimension of the final image to be produced", 512);
        *Ny = input_int("Image height (pixels)", "IMny",
                        "The range dimension in pixels of the final image to be produced", 512);
        *xstart = 0;
        *ystart = 0;
        
        
        *surfaceFilename = input_string("Surface filename", "OutputSurfaceFilename",
                                        "Generated surface will be saved to this file.", "/local_storage/DGM/surface.dat");
        
        // Generate a flat focal surface and save it to file. This will be read back in later
        //
        create_surface(sc_pos, hdr.pulses[*startPulse + *nPulses/2].sat_ps_tx, *surfaceFilename, *Nx, *Ny, surface_sx, surface_sy, status);
        
    }else{
        // Just load up an existing surface - and assume the user knows what they're doing
        // (ie. its of the right area, has sensible pixel size, etc etc).
        //
        *surfaceFilename = input_string("Surface filename", "InputSurfaceFilename",
                                        "The filename of the surface to read in (in Dstl .dat format)",
                                        "/local_storage/DGM/surface.dat");
        SPImage surface ;
        
        im_init(&surface, status);
        im_load_metadata(&surface, *surfaceFilename, status);
        CHECK_STATUS(status) ;
        
        *Nx = (int)surface.nx ;
        *Ny = (int)surface.ny ;
        
        *Nx = input_int("Number of surface samples to process in x", "Nx",
                        "Number of cross-range samples in final image", *Nx) ;
        *Ny = input_int("Number of surface samples to process in y", "Ny",
                        "Number of range samples in final image", *Ny) ;
        
        *xstart = (int)(surface.nx - *Nx)/2 ;
        *ystart = (int)(surface.ny - *Ny)/2 ;
        *xstart = input_int("First surface sample in x for SAR image", "xstart",
                            "Index of first sample in x ditrection within the surface for the first cross-range pixel in the SAR image",
                            *xstart);
        *ystart = input_int("First surface sample in y for SAR image", "ystart",
                            "Index of first sample in y ditrection within the surface for the first range pixel in the SAR image",
                            *ystart);

        getPulses(&hdr, nPulses, startPulse);

    }
    
    int PerformRCSCorr = 0;
    //PerformRCSCorr = input_yesno("Perform RCS correction on image?", "RCSCorr",
    //                            "If yes then each pixel in the output image contains the Radar Cross Section at that location specified in metres^2", PerformRCSCorr);
    if( PerformRCSCorr){
        *EIRP = input_dbl("EIRP", "EIRP",
                          "This is the Effective Isotropic Radiated Power in watts. It is required for the RCS \n \
                          calculation. If Sarcastic(tm) has been used to generate the cphd file then it can be found \n \
                          on the output line", *EIRP);
    }else{
        *EIRP = 0.0;
    }
    
    *outputFilename = input_string("Filename for image", "OutputFName",
                                   "Filename for final SAR image formed by TDP",
                                   "/local_storage/DGM/image1.dat");

    return ;
}

