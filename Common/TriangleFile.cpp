//
//  TriangleFile.cpp
//  sarcastic
//
//  Created by Darren Muff on 21/07/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//

#include "TriangleFile.h"
#include "materialProperties.h"
#include <string.h>

TriangleFile::TriangleFile( std::vector<Triangle> tris ) {
    triangles = tris ;
}

void TriangleFile::WriteFile(std::string fname){
    FILE *fp;
    fp = fopen(fname.c_str(), "w");
    if (fp==NULL) {
        printf("Error : Could not open file %s\n",fname.c_str());
        exit(1);
    }
    
    int ntri = (int)triangles.size() ;
    int vertexBytes = sizeof(double) ;
    int MAXMATNAME = MATBYTES ; // bytes
    
    fwrite(&ntri, sizeof(int), 1, fp);
    fwrite(&vertexBytes, sizeof(int),1,fp);
    fwrite(&MAXMATNAME, sizeof(int), 1, fp);
    
    for (int i=0; i<ntri; i++){
        
        fwrite(&i, sizeof(int), 1, fp);
        
        printf("Triangle %d\n",i);
        double x,y,z, *p ;
        
        Triangle tri = triangles[i];
        x = tri.AA.x;
        y = tri.AA.y;
        z = tri.AA.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.BB.x;
        y = tri.BB.y;
        z = tri.BB.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.CC.x;
        y = tri.CC.y;
        z = tri.CC.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.NN.x;
        y = tri.NN.y;
        z = tri.NN.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.area ;
        fwrite(&x, sizeof(double), 1, fp);
        
        p = tri.globalToLocalMat;
        for (int j=0; j<9; j++){
            x = p[j] ;
            fwrite(&x, sizeof(double), 1, fp);
        }
        p = tri.localToGlobalMat ;
        for (int j=0; j<9; j++){
            x = p[j] ;
            fwrite(&x, sizeof(double), 1, fp);
        }
        char matname[MAXMATNAME];
        strcpy(matname, materialProperties[tri.matId].matname);
        fwrite(matname, sizeof(char), MAXMATNAME, fp);
        
    }
    fclose(fp);


}

//std::ostream & operator<<(std::ostream & os, const TriangleFile & trifile){
//    os.setf(std::ios_base::showpoint);
//
//    os << "Triangle file       : "<<TriangleFile::version<< std::endl;
//    os << "Number of Triangles : "<< TriangleFile::triangles.size();
//    
//}