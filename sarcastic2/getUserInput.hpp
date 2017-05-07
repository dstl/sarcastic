//
//  getUserInput.hpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef getUserInput_hpp
#define getUserInput_hpp

#include <stdio.h>
#include "TriangleMesh.hpp"
#include <SILib2/SIlib2.h>
#include "threadCore.hpp"

int getUserInput(CPHDHeader *hdr, TriangleMesh *baseMesh, TriangleMesh *moverMesh, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam,
                 int *interrogate, SPVector *interogPt, double *interogRad,
                 FILE **interrogateFP, int * pulseUndersampleFactor, SPStatus *status) ;

#endif /* getUserInput_hpp */
