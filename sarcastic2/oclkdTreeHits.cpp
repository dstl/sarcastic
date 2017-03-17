
#include "oclkdTreeHits.hpp"


void oclkdTreeHits(cl_context         context,            // OpenCL context - already instantiated
                   cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                   cl_kernel          STkernel,           // stacklessTraverse kernel, already compiled and created from a program
                   int                nRays,              // Total number of rays to cast
                   size_t             localWorkSize,      // Work dimensions for this device
                   cl_mem             dTriangles,         // device memory for triangles
                   cl_mem             dKdTree,            // device memory for KdTree
                   AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                   Ray *              rays,               // Array of rays (size nRays). Each ray will be cast through KdTree
                   Hit *              hits                // output array of hit locations
)
{
    cl_mem dHits ;
    cl_mem dRays ;
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    // zero hits array so that we can see which rays hit something
    //
    memset(hits, 0, sizeof(Hit)*nRays);
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    dHits        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Hit)*nRays, NULL, &_err));
    dRays        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Ray)*nRays,  NULL, &_err));
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays, CL_TRUE, 0, sizeof(Ray)*nRays, rays, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dHits, CL_TRUE, 0, sizeof(Hit)*nRays, hits, 0, NULL,NULL));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(STkernel, 0,  sizeof(cl_mem), &dKdTree));
    CL_CHECK(clSetKernelArg(STkernel, 3,  sizeof(cl_mem), &dTriangles));
    CL_CHECK(clSetKernelArg(STkernel, 4,  sizeof(AABB),   &SceneBoundingBox));
    CL_CHECK(clSetKernelArg(STkernel, 5,  sizeof(int),    &nRays));
    CL_CHECK(clSetKernelArg(STkernel, 6,  sizeof(cl_mem), &dRays));
    CL_CHECK(clSetKernelArg(STkernel, 7,  sizeof(cl_mem), &dHits));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, STkernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dRays, CL_TRUE,0,sizeof(Ray)*nRays,rays,0,NULL,NULL));
    CL_CHECK(clEnqueueReadBuffer(Q,dHits, CL_TRUE,0,sizeof(Hit)*nRays,hits,0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dHits);
    clReleaseMemObject(dRays);
    
    return ;
}
