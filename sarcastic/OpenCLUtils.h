/** @file********************************************************************
 *
 *       Module:    OpenCLUtils.h
 *      Program:    SIlib
 *   Created by:    Darren Muff on 20/08/2012.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions needed to interpolate images
 *      Where return values are usefu these are passed back as
 *      SPStatus->status ints and the SPStatus parameter set accordingly
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


#ifndef OpenCLUtils_h
#define OpenCLUtils_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#include <SIlib/SIlib.h>
#ifndef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
#define CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV       0x4000
#define CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV       0x4001
#define CL_DEVICE_REGISTERS_PER_BLOCK_NV            0x4002
#define CL_DEVICE_WARP_SIZE_NV                      0x4003
#define CL_DEVICE_GPU_OVERLAP_NV                    0x4004
#define CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV            0x4005
#define CL_DEVICE_INTEGRATED_MEMORY_NV              0x4006
#endif
#define OPENCL_NO_PLATFORM_FOUND        (134)
#define OPENCL_MALLOC_FAIL              (135)
#define OPENCL_NO_DEVICES_FOUND         (136)
#define OPENCL_WORKSIZE_ERROR           (137)
#define OPENCL_X_DIMENSION_ERROR        (138)
#define OPENCL_Y_DIMENSION_ERROR        (139)

/// Error checking wrapper for OpenCl routines that return an error code
///
#define CL_CHECK(_expr)                                                             \
    do {                                                                            \
        cl_int _err = _expr;                                                        \
        if (_err == CL_SUCCESS)                                                     \
            break;                                                                  \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);    \
        abort();                                                                    \
    } while (0)

/// Error checkng wrapper for OpenCl routines that return pass an error parameter
/// into the routine. Make sure the parameter passed into the OpenCL routine
/// is '_err' and is not defined elsewhere
///
#define CL_CHECK_ERR(_expr)                                                         \
({                                                                                  \
    cl_int _err = CL_INVALID_VALUE;                                                 \
    typeof(_expr) _ret = _expr;                                                     \
    if (_err != CL_SUCCESS) {                                                       \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);    \
        abort();                                                                    \
    }                                                                               \
    _ret;                                                                           \
})

typedef struct OCLPlatform {
    int                     nDevs;          ///< Number of opencl devices
    cl_device_id            *device_ids;    ///< array of size 'devsToUse' containing device IDs
    cl_context_properties   *props;         ///< properties for selected platform id
    cl_platform_id clSelectedPlatformID;    ///< Platform id. useful for getting parameters

} OCLPlatform ;

typedef struct OCLTask {
    int                     devsToUse;      ///< Number of opencl devices
    cl_device_id            *device_ids;    ///< array of size 'devsToUse' containing device IDs
    cl_context_properties   *props;         ///< properties for selected platform id
    cl_context              *contexts;      ///< array of size 'devsToUse' containing contexts
    cl_program              *programs;      ///< array of size 'devsToUse' containing programs
    cl_kernel               *kernels;       ///< array of size 'devsToUse' containin kernels
    cl_command_queue        *commands;      ///< array of size 'devsToUse' containing command queues
    const char              *sourceCode;    ///< string wth kernel source code for each device
    const char              *kernelName;    ///< the kernel name (entry point) within the kernel
} OCLTask ;

/// NVIDIA devices have firmware that lets a call know the major and minor
/// compute capabilities of each card.
/// This function allows only cards above a certain major and minor compute capability to be selected
/// This is useful if you don't want to select the card that the monitor is connected to
///
int oclGetNvidiaDevices(OCLPlatform  *platform,             ///< Platform to look for devices
                        int              major,             ///< Minimum NVidia Major compute capability
                        int              minor,             ///< Mimium NVidia Minor compute capability
                        cl_device_id **devices,             ///< Array of devices to be return that meet reqs
                        cl_uint   *deviceCount,             ///< Number of devices that meet reqs
                        SPStatus       *status              ///< SPStatus - error checking and return values
                        );

/// get opencl devices based upon the number of compute units the device has
/// This is useful if you do not want to select the card that is connected to the
/// monitor (and you are not using NVidia cards
///
int oclGetGPUDevices(OCLPlatform          *platform,        ///< Platform to look for devices
                     int min_required_compute_units,        ///< minimum compute units required for a device
                     cl_device_id         **devices,        ///< Array of devices to be return that meet reqs
                     cl_uint           *deviceCount,        ///< Number of devices that meet reqs
                     SPStatus               *status         ///< SPStatus - error checking and return values
                    );

/// Get opencl devices based upon the vendor name and device name
/// if the vendor name is blank then any device will be used
/// if the device name is blank and the device passed the vendor test
/// then any device will be selected.
///
int oclGetNamedGPUDevices(OCLPlatform  *platform,           ///< Platform to look for devices
                          char       *vendorName,           ///< vendor name of device to use (blank for any)
                          char       *deviceName,           ///< device name for devices to use (blank if any)
                          cl_device_id **devices,           ///< Array of devices to be return that meet reqs
                          cl_uint   *deviceCount,           ///< Number of devices that meet reqs
                          SPStatus       *status            ///< SPStatus - error checking and return values
                        );

/// function to return the device ids for CPU devices on the current opencl platform
///
int oclGetCPUDevices(OCLPlatform  *platform,                ///< Platform to look for devices
                     cl_device_id **devices,                ///< Array of devices to be return that meet reqs
                     cl_uint   *deviceCount,                ///< Number of devices that meet reqs
                     SPStatus       *status                 ///< SPStatus - error checking and return values
                    );

/// Get the platform id for later use in OCL functions
///
int oclGetPlatformID(OCLPlatform *platform, SPStatus *status) ;

/// Function to calculate the best 2D work size for a device based upon the 2D requirements for
/// the data (X,Y) and the device capabilities. results are returned in locaXSize and localYSize
///
int best2DWorkSize(cl_kernel Kernel, cl_device_id device, long X, long Y, int *localXSize, int *localYSize, SPStatus *status);

/// Simple function to load in OpenCL kernel source from a file
///
char * loadProgramSource(const char *filename, size_t *size);

/// function to work out how best to fit a 2D data array onto a GPU device that
/// has deviceMemorySize storage available for the data.
/// X,Y are the sizes of the 2D data that will be copied to the device
/// sizeOfElement is the size in bytes of a single element in the 2D array
/// deviceMemorySize is the size avaialble per device that the 2D data will be fit into
/// nDevices is the number of devices available
/// xPerDevice,yPerDevice is the x and y dims of the data that will be copied onto each
/// device for each iteration through the algorithm
/// cycles is the number cycles that each device will be looped over in order to fit
/// all the data through the GPU's
///
int dataForDevices2D(int X, int Y, size_t sizeOfElement, cl_long deviceMemorySize, int nDevices, int *xPerDevice, int *yPerDevice, int *cycles, SPStatus *status);

/// Approach is to split the data into slices (or blocks) in X-direction. Each
/// block has the same number of X-elements as X. If X*sizeOfElement > deviceMemorySize
/// then an error is returned. yPerDevice is the total Y-dimension divided by the number of devices.
/// The number of rows for each block is calculated as
/// the largest number of rows that will fit into deviceMemorySize. This is returned
/// as yPerBlock. If yPerDevice % yPerBlock != 0 then the the remainder is returned in
/// remainder. Cycles is the number of blocks (inc remainder) required for all Y
///
int dataForDevices(int X, int Y, size_t sizeOfElement, cl_long deviceMemorySize, int nDevices, int *yPerdevice, int *yPerBlock, int *remainder, int *cycles, SPStatus *status) ;

/// Take an OpenCl platform definition, a number of devices to use,
/// and a kernel name and build a task structure that contains information
/// about contexts, a loaded kernel and a command queue
///
int oclCreateTask(OCLPlatform *platform, int devsToUse, char *kernelName, const char *kernelCode, OCLTask *task, SPStatus *status );

/// Destroy an OpenCl task structure previously created with oclCreateTask.
/// The function frees up all previously allocated memory.
///
void oclDestroyTask(OCLTask *task, SPStatus *status);
    
void oclPrintPlatform(OCLPlatform platform);

#endif
