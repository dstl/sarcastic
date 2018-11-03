/***************************************************************************
 * 
 *           Module :  sceneExtent.cpp
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 28/05/2018
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :
 *      function to calculate the maximum and minimum extent of a scene when viewed
 *      from a specific location
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

#include "sceneExtent.hpp"
#include <sarclib/sarclib.h>
#include "TriangleMesh.hpp"

void sceneExtent(SPVector Pos, TriangleMesh &mesh, double &maxEl, double &maxAz, double &minEl, double &minAz, AABB &SceneBoundingBox){
    double centreRange ;
    SPVector rVect, unitBeamAz, unitBeamEl, zHat;
    SPVector min,max; VECT_CREATE(9e10, 9e10, 9e10, min); VECT_CREATE(-9e10, -9e10, -9e10, max);
    SPVector boxPts[8] ;
    
    VECT_CREATE(0., 0., 1., zHat);
    centreRange = VECT_MAG(Pos);
    VECT_MINUS( Pos, rVect ) ;
    VECT_CREATE(0, 0, 1., zHat) ;
    VECT_CROSS(rVect, zHat, unitBeamAz);
    VECT_NORM(unitBeamAz, unitBeamAz) ;
    VECT_CROSS(unitBeamAz, rVect, unitBeamEl) ;
    VECT_NORM(unitBeamEl, unitBeamEl) ;
    maxEl = maxAz = minEl = minAz = 0.0 ;
    for(int i=0; i<mesh.triangles.size(); ++i){
        for(int j=0; j<3; ++j){
            min.cell[j] = (mesh.AABBs[i].AA.cell[j] < min.cell[j]) ? mesh.AABBs[i].AA.cell[j] : min.cell[j] ;
            max.cell[j] = (mesh.AABBs[i].BB.cell[j] > max.cell[j]) ? mesh.AABBs[i].BB.cell[j] : max.cell[j] ;
        }
    }
    SceneBoundingBox.AA = min ; SceneBoundingBox.BB = max ;
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[0]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[1]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[2]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[3]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[4]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[5]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[6]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[7]);

    for( int k=0; k<8; k++){
        double El = VECT_DOT(boxPts[k], unitBeamEl) ;
        double Az = VECT_DOT(boxPts[k], unitBeamAz) ;
        maxEl = ( maxEl < El ) ? El : maxEl ;
        maxAz = ( maxAz < Az ) ? Az : maxAz ;
        minEl = ( minEl > El ) ? El : minEl ;
        minAz = ( minAz > Az ) ? Az : minAz ;
    }
    
    return ;
}

