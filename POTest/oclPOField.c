//
//  oclPOField.c
//  sarcastic
//
//  Created by Darren Muff on 08/12/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib/SIlib.h>
#include "OpenCLUtils.h"
#include "POTriangle.h"

int buildKernel(cl_context context,             // OpenCL Context. Already created
                const char *kernelCodePath,     // String defining FQPN for kernel code
                const char *kernelCodeName,     // OpenCL kernel name for kernel code
                cl_device_id devID,             // Device Id on this platform to compile for
                int        workDims,            // WorkDimensions to optimise localWorkSize for
                cl_program *program,            // Output - Compiled OCL programme
                cl_kernel  *kernel,             // Output - OpenCL Kernel
                size_t     *localWorkSize       // Output - Prefered local worksize for kernel
);

void oclPOField(triangle *tris, int ntris, Ray *rays, int nrays, SPVector *hitPoints, SPVector RxPnt, double k, SPVector RXVdir, SPVector RXHdir, SPCmplx *EsV, SPCmplx *EsH)
{
    
    cl_context context ;
    cl_command_queue commandQ ;
    OCLPlatform platform;
    cl_int err ;
    SPStatus status ;
    im_init_status(status, 0);
    
    
    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
//    oclPrintPlatform(platform);
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
    
    // Find number of devices on this platform
    //
    int useGPU = 0 ;
    cl_uint ndevs ;
    
    if (useGPU){
        err = oclGetNamedGPUDevices(&platform, "", "",&platform.device_ids , &ndevs, &status);
        if(err != CL_SUCCESS){
            char cbuf[1024];
            cl_uint max_compute_units;
            printf("GPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE VENDOR             : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
                printf("  DEVICE MAX COMPUTE UNITS  : %d\n",max_compute_units);
            }
        }else{
            printf("No GPU devices available\n");
            exit(1);
        }
    }else{
        err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);
        if(err != CL_SUCCESS){
            printf("CPU DEVICES                 : %d\n",ndevs);
        }
    }
    platform.nDevs = ndevs ;

    // Create OpenCL context and command queue for this device
    //
    cl_device_id devId ;
    devId = platform.device_ids[0];
    
    context  = CL_CHECK_ERR(clCreateContext(platform.props, 1, &devId, NULL, NULL, &_err));
    commandQ = CL_CHECK_ERR(clCreateCommandQueue(context, devId, 0, &_err));
    
    cl_program POFieldPG ;
    cl_kernel  POFieldKL ;
    size_t     POFieldLWS ;
    
    static char *POFieldCode      = "/Users/darren/Development/sarcastic/sarcastic/kernels/POField.cl" ;
    
    CL_CHECK(buildKernel(context, POFieldCode, "POField", devId, 1, &POFieldPG, &POFieldKL, &POFieldLWS));

    
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nrays ;
    while (globalWorkSize % POFieldLWS)globalWorkSize++;
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    cl_mem dtris;
    cl_mem drays;
    cl_mem dhits;
    cl_mem dEsVs;
    cl_mem dEsHs;
    
    dtris = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(triangle)*ntris, NULL, &_err));
    drays = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Ray)*nrays, NULL, &_err));
    dhits = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(SPVector)*nrays, NULL, &_err));
    dEsVs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SPCmplx)*nrays, NULL, &_err));
    dEsHs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SPCmplx)*nrays, NULL, &_err));
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtris, CL_TRUE, 0, sizeof(triangle)*ntris, tris,      0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, drays, CL_TRUE, 0, sizeof(Ray)*nrays,      rays,      0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dhits, CL_TRUE, 0, sizeof(SPVector)*nrays, hitPoints, 0, NULL, NULL));

    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(POFieldKL, 0,   sizeof(cl_mem),   &dtris));
    CL_CHECK(clSetKernelArg(POFieldKL, 1,   sizeof(int),      &ntris));
    CL_CHECK(clSetKernelArg(POFieldKL, 2,   sizeof(cl_mem),   &drays));
    CL_CHECK(clSetKernelArg(POFieldKL, 3,   sizeof(int),      &nrays));
    CL_CHECK(clSetKernelArg(POFieldKL, 4,   sizeof(cl_mem),   &dhits));
    CL_CHECK(clSetKernelArg(POFieldKL, 5,   sizeof(SPVector), &RxPnt));
    CL_CHECK(clSetKernelArg(POFieldKL, 6,   sizeof(double),   &k));
    CL_CHECK(clSetKernelArg(POFieldKL, 7,   sizeof(SPVector), &RXVdir));
    CL_CHECK(clSetKernelArg(POFieldKL, 8,   sizeof(SPVector), &RXHdir));
    CL_CHECK(clSetKernelArg(POFieldKL, 9,   sizeof(cl_mem),   &dEsVs));
    CL_CHECK(clSetKernelArg(POFieldKL, 10,  sizeof(cl_mem),   &dEsHs));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(commandQ, POFieldKL, 1, NULL, &globalWorkSize, &POFieldLWS, 0, NULL, NULL));
    
    // Read results from device
    //
    SPCmplx *VsTmp, *HsTmp;
    VsTmp = sp_malloc(sizeof(SPCmplx) * nrays);
    HsTmp = sp_malloc(sizeof(SPCmplx) * nrays);
    
    CL_CHECK(clEnqueueReadBuffer(commandQ,dEsVs, CL_TRUE,0,sizeof(SPCmplx)*nrays, VsTmp, 0,NULL,NULL));
    CL_CHECK(clEnqueueReadBuffer(commandQ,dEsHs, CL_TRUE,0,sizeof(SPCmplx)*nrays, HsTmp, 0,NULL,NULL));
    
    SPCmplx VsSum, HsSum ;
    VsSum.r = VsSum.i = HsSum.r = HsSum.i = 0.0f ;
    
    for (int r=0; r<nrays; r++ ){
        
        CMPLX_ADD(VsTmp[r], VsSum, VsSum) ;
        CMPLX_ADD(HsTmp[r], HsSum, HsSum) ;
    }
    
    *EsV = VsSum ;
    *EsH = HsSum ;
    
    free(VsTmp) ;
    free(HsTmp) ;
    
    clReleaseMemObject(dtris);
    clReleaseMemObject(drays);
    clReleaseMemObject(dhits);
    clReleaseMemObject(dEsVs);
    clReleaseMemObject(dEsHs);
    clReleaseKernel(POFieldKL) ;
    clReleaseCommandQueue(commandQ) ;
    clReleaseContext(context) ;
    clReleaseProgram(POFieldPG) ;
    
    return ;
}
