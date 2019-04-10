/***************************************************************************
 * 
 *           Module :  OpenCLUtils.c
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 20/08/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      This file contains useful finctions for setting up OpenCL platforms
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

#include "OpenCLUtils.h"

// NVIDIA devices have firmware that lets a call know the major and minor
// compute capabilities of each card.
// This function allows only cards above a certain major and minor compute capability to be selected
// This is useful if you don't want to select the card that the monitor is connected to
//
int oclGetNvidiaDevices(OCLPlatform *platform,      // Platform ID to look for devices
                        int major,                  // Minimum NVidia Major compute capability
                        int minor,                  // Mimium NVidia Minor compute capability
                        cl_device_id **devices,     // Array of devices to be return that meet reqs
                        cl_uint *deviceCount,       // Number of devices that meet reqs
                        SPStatus *status            // SPStatus - error checking and return values
                        ){
    cl_uint ciErrNum;
    CHECK_STATUS(status);
    ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, 0, NULL, deviceCount);
    // check for 0 devices found or errors...
    //
    if (deviceCount == 0){
        if(status->debug >= 10 )printf(" No devices found supporting OpenCL (return code %i)\n\n", ciErrNum);
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return( status->status );
    }else if (ciErrNum != CL_SUCCESS){
        if(status->debug >= 10 )printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
        status->status = ciErrNum ;
        return( status->status );
    }else{
        if ((*devices = (cl_device_id*)malloc(sizeof(cl_device_id) * *deviceCount)) == NULL){
            if(status->debug >= 10 )printf(" Failed to allocate memory for devices !!!\n\n");
            status->status = OPENCL_MALLOC_FAIL ;
            return( status->status );
        }
        ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, *deviceCount, *devices, NULL);
        if(ciErrNum != CL_SUCCESS){
            if(status->debug >= 10 )printf(" Failed to collect GPU device IDs. Error %i\n",ciErrNum);
            status->status = ciErrNum ;
            return( status->status );
        }
        int selectdevcount =0;
        for(int i=0; i< *deviceCount; i++){
            cl_uint compute_capability_major, compute_capability_minor;
            ciErrNum = clGetDeviceInfo((*devices)[i], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &compute_capability_major, NULL);
            ciErrNum |= clGetDeviceInfo((*devices)[i], CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, sizeof(cl_uint), &compute_capability_minor, NULL);
            if(ciErrNum != CL_SUCCESS){
                if(status->debug >= 10 )printf(" Failed to collect device IDs for COMPUTE_CAPABILITY. Error %i\n",ciErrNum);
                status->status = ciErrNum ;
                return( status->status );
            }
            if(compute_capability_major >= major && compute_capability_minor >= minor){
                (*devices)[selectdevcount] = (*devices)[i];
                selectdevcount++;
            }
        }
        *deviceCount         = selectdevcount;
        platform->nDevs      = selectdevcount ;
        platform->device_ids = *devices ;
    }
    
    return( status->status );
}

// get opencl devices based upon the number of compute units the device has
// This is useful if you do not want to select the card that is connected to the
// monitor (and you are not using NVidia cards
//
int oclGetGPUDevices(OCLPlatform          *platform,    // Platform ID to look for devices
                     int min_required_compute_units,    // minimum compute units required for a device
                     cl_device_id         **devices,    // Array of devices to be return that meet reqs
                     cl_uint           *deviceCount,    // Number of devices that meet reqs
                     SPStatus               *status     // SPStatus - error checking and return values
                    ){
    cl_uint ciErrNum;
    CHECK_STATUS(status) ;

    ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, 0, NULL, deviceCount);
    // check for 0 devices found or errors...
    if (deviceCount == 0){
        if(status->debug >= 10)printf(" No GPU devices found supporting OpenCL (return code %i)\n\n", ciErrNum);
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return( status->status );
    }else if (ciErrNum != CL_SUCCESS){
        if(status->debug >= 10)printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
        status->status = ciErrNum ;
        return( status->status );
    }else{
        if ((*devices = (cl_device_id*)malloc(sizeof(cl_device_id) * *deviceCount)) == NULL){
            if(status->debug >= 10)printf(" Failed to allocate memory for devices !!!\n\n");
            status->status = OPENCL_MALLOC_FAIL ;
            return( status->status );
        }
        ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, *deviceCount, *devices, NULL);
        if(ciErrNum != CL_SUCCESS){
            if(status->debug >= 10)printf(" Failed to collect GPU device IDs. Error %i\n",ciErrNum);
            status->status = ciErrNum ;
            return( status->status );
        }
        int selectdevcount =0;
        for(int i=0; i< *deviceCount; i++){
            cl_uint max_compute_units;
            ciErrNum = clGetDeviceInfo((*devices)[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
            if(ciErrNum != CL_SUCCESS){
                if(status->debug >= 10)printf(" Failed to collect device IDs for COMPUTE_CAPABILITY. Error %i\n",ciErrNum);
                status->status = ciErrNum ;
                return( status->status );
            }
            if(max_compute_units >= min_required_compute_units){
                (*devices)[selectdevcount] = (*devices)[i];
                selectdevcount++;
            }
        }
        *deviceCount         = selectdevcount;
        platform->nDevs      = selectdevcount ;
        platform->device_ids = *devices ;
    }
    
    return( status->status );
    
}

// Get opencl devices based upon the vendor name and device name
// if the vendor name is blank then any device will be used
// if the device name is blank and the device passed the vendor test
// then any device will be selected.
//
int oclGetNamedGPUDevices(OCLPlatform  *platform,       // Platform to look for devices
                          char       *vendorName,       // vendor name of device to use (blank for any)
                          char       *deviceName,       // device name for devices to use (blank if any)
                          cl_device_id **devices,       // Array of devices to be return that meet reqs
                          cl_uint   *deviceCount,       // Number of devices that meet reqs
                          SPStatus       *status        // SPStatus - error checking and return values
                        ){
    cl_uint ciErrNum;
    CHECK_STATUS(status) ;
    
    ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, 0, NULL, deviceCount);
    // check for 0 devices found or errors...
    if (deviceCount == 0){
        if(status->debug >= 10)printf(" No GPU devices found supporting OpenCL (return code %i)\n\n", ciErrNum);
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return( status->status );
    }else if (ciErrNum != CL_SUCCESS){
        if(status->debug >= 10)printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
        status->status = ciErrNum ;
        return( status->status );
    }else{
        if ((*devices = (cl_device_id*)malloc(sizeof(cl_device_id) * *deviceCount)) == NULL){
            if(status->debug >= 10)printf(" Failed to allocate memory for devices !!!\n\n");
            status->status = OPENCL_MALLOC_FAIL ;
            return( status->status );
        }
        ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_GPU, *deviceCount, *devices, NULL);
        if(ciErrNum != CL_SUCCESS){
            if(status->debug >= 10)printf(" Failed to collect GPU device IDs. Error %i\n",ciErrNum);
            status->status = ciErrNum ;
            return( status->status );
        }
        int selectdevcount =0;
        for(int i=0; i< *deviceCount; i++){
            char cbuf[1024];
            ciErrNum = clGetDeviceInfo((*devices)[i], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
            if(strlen(vendorName) == 0 || !strcmp(cbuf, vendorName)){
                ciErrNum |= clGetDeviceInfo((*devices)[i], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                if (strlen(deviceName) == 0 || !strcmp(cbuf, deviceName)) {
                    (*devices)[selectdevcount] = (*devices)[i];
                    selectdevcount++;
                }
                
            }
            if(ciErrNum != CL_SUCCESS){
                if(status->debug >= 10)printf(" Failed to collect device IDs for COMPUTE_CAPABILITY. Error %i\n",ciErrNum);
                status->status = ciErrNum ;
                return( status->status );
            }
        }
        *deviceCount         = selectdevcount;
        platform->nDevs      = selectdevcount ;
        platform->device_ids = *devices ;
    }
    
    return( status->status );
}

// function to return the device ids for CPU devices on the current opencl platform
//
int oclGetCPUDevices(OCLPlatform  *platform,    // Platform to look for devices
                     cl_device_id **devices,    // Array of devices to be return that meet reqs
                     cl_uint   *deviceCount,    // Number of devices that meet reqs
                     SPStatus       *status     // SPStatus - error checking and return values
                    ){
    cl_uint ciErrNum;
    CHECK_STATUS(status) ;

    ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_DEFAULT, 0, NULL, deviceCount);
    // check for 0 devices found or errors...
    //
    if (deviceCount == 0){
        if(status->debug >= 10)printf(" No devices found supporting OpenCL (return code %i)\n\n", ciErrNum);
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return( status->status );
    }else if (ciErrNum != CL_SUCCESS){
        printf("Device count %d\n", *deviceCount);
        if(status->debug >= 10)printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
        status->status = ciErrNum ;
        return( status->status );
    }else{
        if ((*devices = (cl_device_id*)malloc(sizeof(cl_device_id) * *deviceCount)) == NULL){
            if(status->debug >= 10)printf(" Failed to allocate memory for devices !!!\n\n");
            status->status = OPENCL_MALLOC_FAIL ;
            return( status->status );
        }
        ciErrNum = clGetDeviceIDs (platform->clSelectedPlatformID, CL_DEVICE_TYPE_DEFAULT, *deviceCount, *devices, NULL);
        if(ciErrNum != CL_SUCCESS){
            if(status->debug >= 10)printf(" Failed to collect GPU device IDs. Error %i\n",ciErrNum);
            status->status = ciErrNum ;
            return( status->status );
        }
        platform->nDevs = *deviceCount ;
    }
    
    return( status->status );
}

// Get the platform id for later use in OCL functions
//
int oclGetPlatformID(OCLPlatform *platform, SPStatus *status)
{
    char chBuffer[1024];
    cl_uint num_platforms;
    cl_platform_id* clPlatformIDs;
    cl_int ciErrNum;
    platform->clSelectedPlatformID = NULL;
    
    CHECK_STATUS(status) ;
    
    // Get OpenCL platform count
    //
    ciErrNum = clGetPlatformIDs (0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS){
        if(status->debug > 10)printf(" Error %i in clGetPlatformIDs Call !!!\n\n", ciErrNum);
        status->status=ciErrNum ;
        return status->status ;
    }else{
        if(num_platforms == 0){
            if(status->debug > 10)printf("No OpenCL platform found!\n\n");
            status->status = OPENCL_NO_PLATFORM_FOUND ;
            return status->status ;
        }else{
            // if there's a platform or more, make space for ID's
            //
            if ((clPlatformIDs = (cl_platform_id*)malloc(num_platforms * sizeof(cl_platform_id))) == NULL){
                if(status->debug > 10)printf("Failed to allocate memory for cl_platform ID's!\n\n");
                status->status = OPENCL_MALLOC_FAIL ;
                return status->status ;
            }
            
            // get platform info for each platform and trap the NVIDIA platform if found
            //
            ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);
            if (ciErrNum != CL_SUCCESS) {
                if(status->debug > 10)printf("ERROR: Failed to get platform IDs\n");
                status->status = ciErrNum ;
                free(clPlatformIDs);
                return status->status ;
            }
            for(cl_uint i = 0; i < num_platforms; ++i){
                ciErrNum = clGetPlatformInfo (clPlatformIDs[i], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);
                if(ciErrNum == CL_SUCCESS){
                    if(strstr(chBuffer, "NVIDIA") != NULL){
                        platform->clSelectedPlatformID = clPlatformIDs[i];
                        break;
                    }
                }
            }
            
            // default to zeroeth platform if NVIDIA not found
            //
            if(platform->clSelectedPlatformID == NULL)
            {
                if(status->debug > 10)printf("WARNING: NVIDIA OpenCL platform not found - defaulting to first platform!\n");
                platform->clSelectedPlatformID = clPlatformIDs[0];
            }
            
            free(clPlatformIDs);
        }
    }
    return CL_SUCCESS;
}

// Function to calculate the best 2D work size for a device based upon the 2D requirements for
// the data (X,Y) and the device capabilities. results are returned in locaXSize and localYSize
//
int best2DWorkSize(cl_kernel Kernel, cl_device_id device, long X, long Y, int *localXSize, int *localYSize, SPStatus *status){
    size_t localWorkSize;
    cl_int err;
    cl_uint maxDims;
    size_t *maxSizes;
    CHECK_STATUS(status);
    
    err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &maxDims, NULL);
    if (err != CL_SUCCESS){
        if(status->debug>=10)printf("Error: Failed to retrieve kernel work group Size! %d\n", err);
        status->status = err ;
        return ( status->status ) ;
    }
    
    maxSizes = (size_t *)malloc(sizeof(size_t)*maxDims);
    if(maxSizes == NULL){
        if(status->debug>=10)printf("ERROR : Malloc failed\n");
        status->status = OPENCL_MALLOC_FAIL ;
        return ( status->status ) ;
    }
    
    err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*maxDims, maxSizes, NULL);
    if (err != CL_SUCCESS){
        if(status->debug>=10)printf("Error: Failed to retrieve kernel work group Size! %d\n", err);
        status->status = err ;
        free(maxSizes);
        return ( status->status );
    }
    
    if ( status->debug >=10 ){
        printf("Max work items per dimension : %zd",maxSizes[0]);
        for(int i=1; i<maxDims; i++)printf(" / %zd",maxSizes[i]);
        printf("\n");
    }
    
    err = clGetKernelWorkGroupInfo(Kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(localWorkSize), &localWorkSize, NULL);
    if (err != CL_SUCCESS){
        if(status->debug>=10)printf("Error: Failed to retrieve kernel work group info! %d\n", err);
        status->status = err ;
        free(maxSizes);
        return ( status->status );
    }

    size_t warp;
    err = clGetKernelWorkGroupInfo(Kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(warp), &warp, NULL);
    if (err != CL_SUCCESS){
        if(status->debug>=10)printf("Error: Failed to retrieve kernel warp info! %d\n", err);
        status->status = err ;
        free(maxSizes);
        return ( status->status );
    }
    
    localWorkSize = (size_t)sqrtf(localWorkSize);
    // Find the maximum localSize that verifies the following conditions:
    // - localXSize * localYSize <= CL_KERNEL_WORK_GROUP_SIZE (hardware limitation)
    // - localXSize,localYSize divides the globalWorkSize's width and height respectively (OpenCL requirement)
    // - device.alignment() divides localXSize * localYSize in order to maximize
    //   processor usage. NVIDIA calls this defaultAligment "warp size".
    // The final value is usually 8 or 16
    //
    *localXSize = (int)localWorkSize;
    
    while ( (*localXSize > 0) && ((X % *localXSize) || (Y % *localXSize) || ((*localXSize * *localXSize) % warp))) {
        
        (*localXSize)--;
    }
    
    if (! *localXSize) {
        // Couldn't find a proper localSize, chances are that globalWorkSize's sides are not
        // divisible by 16, 8, etc.
        // Ignore alignment and try again, in the worst case localSide will be 1
        //
        *localXSize= (int)localWorkSize;
        while ((X % *localXSize) || (Y % *localXSize))
            (*localXSize)--;
    }
    
    if ( *localXSize > maxSizes[0]){
        if(status->debug>=10)printf("Error : max work items per dimension smaller than calculated local X Size\n");
        status->status = OPENCL_WORKSIZE_ERROR ;
        free(maxSizes);
        return ( status->status );
    }
    
    *localYSize = ( maxSizes[1] < *localXSize ) ? (int)maxSizes[1] : *localXSize ;
    
    free(maxSizes);
    return ( status->status );
}

// Simple function to load in OpenCL kernel source from a file
//
char * loadProgramSource(const char *filename, size_t *size){
	
	struct stat statbuf;
	FILE *fh;
	char *source;
	
	fh = fopen(filename, "r");
	if (fh == 0){
        printf("Error : Failed to load CL code from file %s\n",filename);
        exit(0);
    }
	
	stat(filename, &statbuf);
	source = (char *) malloc(statbuf.st_size + 1);
	fread(source, statbuf.st_size, 1, fh);
	source[statbuf.st_size] = '\0';
	*size = statbuf.st_size + 1;
    
    fclose(fh);
	return source;
}

// function to work out how best to fit a 2D data array onto a GPU device that
// has deviceMemorySize storage available for the data.
// X,Y are the sizes of the 2D data that will be copied to the device
// sizeOfElement is the size in bytes of a single element in the 2D array
// deviceMemorySize is the size avaialble per device that the 2D data will be fit into
// nDevices is the number of devices available
// xPerDevice,yPerDevice is the x and y dims of the data that will be copied onto each
// device for each iteration through the algorithm
// cycles is the number cycles that each device will be looped over in order to fit
// all the data through the GPU's
//
int dataForDevices2D(int X, int Y, size_t sizeOfElement,
                      cl_long deviceMemorySize,
                      int nDevices,
                      int *xPerDevice, int *yPerDevice, int *cycles, SPStatus *status)

{
    // Approach is to x-dim of data onto devices first and then work out
    // how many Y lines can be fit. Then loop through all devices until
    // all data has been processed.
    //
    
    if (X % nDevices != 0) {
        if(status->debug>=10)printf("ERROR : X dimension of data is not divisible by %d devices\n",nDevices);
        status->status = OPENCL_X_DIMENSION_ERROR ;
        return ( status->status ) ;
    }
    
    cl_long totalMemAvail = deviceMemorySize * nDevices ;
    
    *xPerDevice = X / nDevices ;
    int xMemSize = *xPerDevice * nDevices * (int)sizeOfElement ;
    
    if(xMemSize > totalMemAvail ){
        if(status->debug>=10)printf("ERROR : one line of 2D data is too large for all devices\n");
        status->status = OPENCL_X_DIMENSION_ERROR ;
        return ( status->status ) ;
    }
    
    *yPerDevice = (int)(totalMemAvail / xMemSize) ;
    *yPerDevice = (*yPerDevice > Y) ? Y : *yPerDevice ;
    *cycles     = Y / *yPerDevice ;
    
    return ( status->status ) ;
    
}


// Approach is to split the data into slices (or blocks) in X-direction. Each
// block has the same number of X-elements as X. If X*sizeOfElement > deviceMemorySize
// then an error is returned. yPerDevice is the total Y-dimension divided by the number of devices.
// The number of rows for each block is calculated as
// the largest number of rows that will fit into deviceMemorySize. This is returned
// as yPerBlock. If yPerDevice % yPerBlock != 0 then the the remainder is returned in
// remainder. Cycles is the number of blocks (inc remainder) required for all Y
// 
int dataForDevices(int X, int Y, size_t sizeOfElement, // number of data elements and size in bytes of each element
                   cl_long deviceMemorySize,                // device memory size in bytes
                   int nDevices,                            // number of OpenCl devices
                   int *yPerDevice,                         // answer. number of rows per device
                   int *yPerBlock,                          // answer. y-size per device.
                   int *remainder,                          // answer. remainder rows per device
                   int *cycles,                             // number of cycles per device
                   SPStatus *status)                         // SPStatus - error checking and return values
{
    if ( X % nDevices != 0) {
        if(status->debug>=10)printf("ERROR : X dimension of data is not divisible by %d devices\n",nDevices);
        status->status = OPENCL_X_DIMENSION_ERROR ;
        return ( status->status ) ;
    }
    if( Y % nDevices != 0 ) {
        if(status->debug>=10)printf("ERROR : Y dimension of data is not divisible by %d devices\n",nDevices);
        status->status = OPENCL_Y_DIMENSION_ERROR ;
        return ( status->status ) ;
    }
    *yPerDevice       = Y / nDevices ;
    cl_long memPerRow = X * sizeOfElement ;
    *yPerBlock        = (int)(deviceMemorySize / memPerRow) ;
    if (*yPerBlock > *yPerDevice ) *yPerBlock = *yPerDevice ;
    *remainder        = *yPerDevice % *yPerBlock ;
    *cycles           = *yPerDevice / *yPerBlock ;
    if (*remainder > 0 ) (*cycles)++;
    
    return ( status->status );
}

// Take an OpenCl platform definition, a number of devices to use,
// and a kernel name and build a task structure that contains information
// about contexts, a loded kernel and a command queue
//
int oclCreateTask(OCLPlatform *platform, int devsToUse, char *kernelName, const char *kernelCode, OCLTask *task, SPStatus *status ){
    
    cl_int err;
    CHECK_STATUS(status);
    
    if(platform->nDevs == 0){
        if(status->debug>=10)printf("ERROR : Platform does not have any available devices ! \n");
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return ( status->status );
    }
    
    task->devsToUse = (devsToUse < platform->nDevs) ? devsToUse : platform->nDevs ;
    if(devsToUse == 0){
        if(status->debug>=10)printf("ERROR : Incorrect number of devices requested !\n");
        status->status = OPENCL_NO_DEVICES_FOUND ;
        return ( status->status ) ;
    }
   
    task->device_ids = (cl_device_id *)sp_malloc(sizeof(cl_device_id)*task->devsToUse);
    task->props      = platform->props ;
    task->contexts   = (cl_context *)sp_malloc(sizeof(cl_context)*(task->devsToUse));
    task->programs   = (cl_program *)sp_malloc(sizeof(cl_program)*(task->devsToUse));
    task->kernels    = (cl_kernel *)sp_malloc(sizeof(cl_kernel)*(task->devsToUse));
    task->commands   = (cl_command_queue *)sp_malloc(sizeof(cl_command_queue)*(task->devsToUse));
    task->kernelName = (char *)sp_calloc(sizeof(char), strlen(kernelName)+1);
    strcpy((char *)task->kernelName, kernelName) ;
    
    for (int dev=0; dev<task->devsToUse; dev++){
        
        task->device_ids[dev] = platform->device_ids[dev] ;
        
        task->contexts[dev] = clCreateContext(task->props, 1, &(task->device_ids[dev]), NULL, NULL, &err);
        if (!task->contexts[dev]){
            if(status->debug>=10)printf("Error [%d] : Failed to create a compute context! \n",err);
            status->status = err ;
            return ( status->status ) ;
        }
        FILE *fp;
        size_t codeSize = 0;
        char *code;
        
        if ( (fp = fopen(kernelCode, "r"))){
            fclose(fp);
            code = loadProgramSource(kernelCode, &codeSize);
            task->programs[dev] = clCreateProgramWithSource(task->contexts[dev], 1, (const char **) &code, &codeSize, &err);
        }else{
            task->programs[dev] = clCreateProgramWithSource(task->contexts[dev], 1, (const char **) &kernelCode, NULL, &err);
        }
        
        if (!task->programs[dev]){
            if(status->debug>=10)printf("Error [%d]: Failed to create compute program for device %d!\n",err,dev);
            status->status = err ;
            return ( status->status ) ;
        }
        err = clBuildProgram(task->programs[dev], 0, NULL, "-Werror", NULL, NULL);
        if (err != CL_SUCCESS){
            size_t len;
            char buffer[32768];
            if(status->debug>=10){
                printf("Error: Failed to build program executable!\n");
                clGetProgramBuildInfo(task->programs[dev], task->device_ids[dev], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
                if(status->debug>=10)printf("err: %d. Buffer:\n",err);
                if(status->debug>=10)printf("%s\n", buffer);
            }
            status->status = err ;
            return ( status->status ) ;
        }
        task->kernels[dev] = clCreateKernel(task->programs[dev], task->kernelName, &err);
        if (!task->kernels[dev] || err != CL_SUCCESS){
            if(status->debug>=10)printf("Error: Failed to build kernel \"%s\" ! : %d\n",task->kernelName,err);
            status->status = err ;
            return ( status->status ) ;
        }
        task->commands[dev] = clCreateCommandQueue(task->contexts[dev], task->device_ids[dev], 0, &err);
        if (!task->commands[dev]){
            if(status->debug>=10)printf("Error: Failed to create a command commands!\n");
            status->status = err ;
            return ( status->status ) ;
        }
    }
    return ( status->status ) ;
}

void oclDestroyTask(OCLTask *task, SPStatus *status){
    CHECK_STATUS(status);
    free(task->device_ids);
    free(task->contexts);
    free(task->programs);
    free(task->kernels);
    free(task->commands);
    free((char *)task->kernelName);
    return ;
}

void oclPrintPlatform(OCLPlatform platform){
    cl_int err;
    char chBuffer[1024];
    
    printf("-----------------------------------------------------------------------------------------\n");
    printf("    OpenCL Device information\n");
    
    err = clGetPlatformInfo (platform.clSelectedPlatformID, CL_PLATFORM_PROFILE, 1024, &chBuffer, NULL);
    if(err == CL_SUCCESS){
        printf("CL_PLATFORM_PROFILE         : %s\n",chBuffer);
    }
    err = clGetPlatformInfo (platform.clSelectedPlatformID, CL_PLATFORM_VERSION, 1024, &chBuffer, NULL);
    if(err == CL_SUCCESS){
        printf("CL_PLATFORM_VERSION         : %s\n",chBuffer);
    }
    err = clGetPlatformInfo (platform.clSelectedPlatformID, CL_PLATFORM_NAME, 1024, &chBuffer, NULL);
    if(err == CL_SUCCESS){
        printf("CL_PLATFORM_NAME            : %s\n",chBuffer);
    }
    err = clGetPlatformInfo (platform.clSelectedPlatformID, CL_PLATFORM_VENDOR, 1024, &chBuffer, NULL);
    if(err == CL_SUCCESS){
        printf("CL_PLATFORM_VENDOR          : %s\n",chBuffer);
    }
    
}








