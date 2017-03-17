//
//  oclkdTreeHits.hpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef oclkdTreeHits_hpp
#define oclkdTreeHits_hpp

#include <stdio.h>
#include "OpenCLUtils.h"
#include "AABB.hpp"
#include "threadCore.hpp"

void oclKdTreeHits(cl_context         context,            // OpenCL context - already instantiated
                   cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                   cl_kernel          STkernel,           // stacklessTraverse kernel, already compiled and created from a program
                   int                nRays,              // Total number of rays to cast
                   size_t             localWorkSize,      // Work dimensions for this device
                   cl_mem             dTriangles,         // device memory for triangles
                   cl_mem             dKdTree,            // device memory for KdTree
                   cl_mem             dtriListData,       // device array containing all indices. size=nTriIndices
                   cl_mem             dtriListPtrs,       // device array containing index in triIndices that matches node. - size=nLeaves
                   AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                   Ray *              rays,               // Array of rays (size nRays). Each ray will be cast through KdTree
                   Hit *              hits                // output array of hit locations
) ;

#endif /* oclkdTreeHits_hpp */
