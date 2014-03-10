//
//  main.c
//  doublePrecisionTest
//
//  Created by Darren on 16/08/2013.
//  Copyright (c) 2013 Dstl. All rights reserved.
//

#include <stdio.h>
#include <OpenCL/opencl.h>
#include <math.h>

//#define NUM_DATA 1

#define CL_CHECK(_expr)                                                         \
    do {                                                                         \
        cl_int _err = _expr;                                                       \
        if (_err == CL_SUCCESS)                                                    \
            break;                                                                   \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);   \
        abort();                                                                   \
} while (0)

#define CL_CHECK_ERR(_expr)                                                     \
({                                                                           \
    cl_int _err = CL_INVALID_VALUE;                                            \
    typeof(_expr) _ret = _expr;                                                \
    if (_err != CL_SUCCESS) {                                                  \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err); \
        abort();                                                                 \
    }                                                                          \
    _ret;                                                                      \
})

void pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
	fprintf(stderr, "OpenCL Error (via pfn_notify): %s\n", errinfo);
}

int main(int argc, const char * argv[])
{
    
    cl_platform_id platforms[100];
	cl_uint platforms_n = 0;
	CL_CHECK(clGetPlatformIDs(100, platforms, &platforms_n));
    
//	printf("=== %d OpenCL platform(s) found: ===\n", platforms_n);
//	for (int i=0; i<platforms_n; i++)
//	{
//		char buffer[10240];
//		printf("  -- %d --\n", i);
//		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, 10240, buffer, NULL));
//		printf("  PROFILE = %s\n", buffer);
//		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 10240, buffer, NULL));
//		printf("  VERSION = %s\n", buffer);
//		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 10240, buffer, NULL));
//		printf("  NAME = %s\n", buffer);
//		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 10240, buffer, NULL));
//		printf("  VENDOR = %s\n", buffer);
//		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, 10240, buffer, NULL));
//		printf("  EXTENSIONS = %s\n", buffer);
//	}
    
	if (platforms_n == 0)
		return 1;
    
	cl_device_id devices[100];
	cl_uint devices_n = 0;
	CL_CHECK(clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_CPU, 100, devices, &devices_n));
    
//	printf("=== %d OpenCL device(s) found on platform:\n", platforms_n);
//	for (int i=0; i<devices_n; i++)
//	{
//		char buffer[10240];
//		cl_uint buf_uint;
//		cl_ulong buf_ulong;
//		printf("  -- %d --\n", i);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL));
//		printf("  DEVICE_NAME = %s\n", buffer);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL));
//		printf("  DEVICE_VENDOR = %s\n", buffer);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL));
//		printf("  DEVICE_VERSION = %s\n", buffer);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL));
//		printf("  DRIVER_VERSION = %s\n", buffer);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL));
//		printf("  DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL));
//		printf("  DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
//		CL_CHECK(clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL));
//		printf("  DEVICE_GLOBAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);
//	}
    
	if (devices_n == 0)
		return 1;
    
	cl_context context;
	context = CL_CHECK_ERR(clCreateContext(NULL, 1, devices, &pfn_notify, NULL, &_err));
    
    double num = M_PI ;
    printf("number before is : %1.17e\n",num);
    
	const char *program_source[] = {
        "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n",
		"__kernel void simple_demo(__global double *src, __global double *dst)\n",
		"{\n",
		"	int i = get_global_id(0);\n",
        "   printf(\"src on device is : %1.17e\\n\",src[i]);\n",
		"	dst[i] = src[i];\n",
		"}\n"
	};
    
	cl_program program;
	program = CL_CHECK_ERR(clCreateProgramWithSource(context, sizeof(program_source)/sizeof(*program_source), program_source, NULL, &_err));
	if (clBuildProgram(program, 1, devices, "", NULL, NULL) != CL_SUCCESS) {
		char buffer[10240];
		clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
		fprintf(stderr, "CL Compilation failed:\n%s", buffer);
		abort();
	}
    
	cl_mem input_buffer;
	input_buffer = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(double), NULL, &_err));
    
	cl_mem output_buffer;
	output_buffer = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double), NULL, &_err));
        
	cl_kernel kernel;
	kernel = CL_CHECK_ERR(clCreateKernel(program, "simple_demo", &_err));
	CL_CHECK(clSetKernelArg(kernel, 0, sizeof(input_buffer), &input_buffer));
	CL_CHECK(clSetKernelArg(kernel, 1, sizeof(output_buffer), &output_buffer));
    
	cl_command_queue queue;
	queue = CL_CHECK_ERR(clCreateCommandQueue(context, devices[0], 0, &_err));
    
    CL_CHECK(clEnqueueWriteBuffer(queue, input_buffer, CL_TRUE, 0, sizeof(double), &num, 0, NULL, NULL));
    
	cl_event kernel_completion;
	size_t global_work_size[1] = { 1 };
	CL_CHECK(clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global_work_size, NULL, 0, NULL, &kernel_completion));
	CL_CHECK(clWaitForEvents(1, &kernel_completion));
	CL_CHECK(clReleaseEvent(kernel_completion));
    
	printf("number after is  :");
    double data;
    CL_CHECK(clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0, sizeof(double), &data, 0, NULL, NULL));
    printf(" %1.17e", data);
	printf("\n");
    
	CL_CHECK(clReleaseMemObject(input_buffer));
	CL_CHECK(clReleaseMemObject(output_buffer));
    
	CL_CHECK(clReleaseKernel(kernel));
	CL_CHECK(clReleaseProgram(program));
	CL_CHECK(clReleaseContext(context));
    
	return 0;
}

