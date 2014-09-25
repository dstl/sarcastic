//
//  main.m
//  colladaToTriFile
//
//  Created by Darren on 21/09/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "MVectorObjc.h"
#import "COLScene.h"
#import "Triangle.h"
#include <SILib.h>

#define ROOTPATH "/Users/Darren/Development"

int main(int argc, const char * argv[])
{
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "colladaToTriFile", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    
    char *str ;
    
    @autoreleasepool {
        
        str = input_string((char *)"Collada filename", (char *)"filename",
                           (char *)"The name of a model in the collada DAE format. Export from Sketchup(R)",
                           (char *) ROOTPATH"/Models/model.dae");
        
        NSString *filestr = [NSString stringWithUTF8String:str];
        if ([filestr rangeOfString:@"file://"].location == NSNotFound) {
            NSMutableString *newFileStr = [NSMutableString stringWithUTF8String:"file://"] ;
            filestr = [newFileStr stringByAppendingString:filestr];
        }
        NSURL * ipurl = [NSURL URLWithString:filestr];
        
        str =input_string((char *)"output filename", (char *)"opFname",
                          (char *)"Full path name for output triangle file",
                          (char *) ROOTPATH"/DATA/triangles.tri");
        
        FILE *fp;
        fp = fopen(str, "w");
        if (fp==NULL) {
            printf("Error : Could not open file %s\n",str);
            exit(1);
        }
        
        /*
         Output File Structure
         int number_of_triangles
         int sizeof_vertex_information
         int length_of_material_names
         foreach triangle{
            int triangle_number
            <sizeof_vertex_information> vertexA.x
            <sizeof_vertex_information> vertexA.y
            <sizeof_vertex_information> vertexA.z
            <sizeof_vertex_information> vertexB.x
            <sizeof_vertex_information> vertexB.y
            <sizeof_vertex_information> vertexB.z
            <sizeof_vertex_information> vertexC.x
            <sizeof_vertex_information> vertexC.y
            <sizeof_vertex_information> vertexC.z
            <sizeof_vertex_information> normal.x
            <sizeof_vertex_information> normal.y
            <sizeof_vertex_information> normal.z
            char[<length_of_material_names>] Material_name
         }
         */
        
        COLScene * scene = [COLScene sceneWithURL:ipurl];
        printf("Scene Summary\n");
        printf("==================\n");
        int ntri = (int)[[scene triangles] count] ;
        fwrite(&ntri, sizeof(int), 1, fp);
        int sizeofdouble = sizeof(double);
        fwrite(&sizeofdouble, sizeof(int),1,fp);
        int nmat = (int)[[scene triangleMaterials] count] ;
        printf("Triangles in scene : %d \n", ntri);
        printf("Materials in scene : %d \n", nmat);
        int MAXMATNAME = 128 ;
        fwrite(&MAXMATNAME, sizeof(int), 1, fp);

        for (int i=0; i<ntri; i++){
           
            fwrite(&i, sizeof(int), 1, fp);
          
            printf("Triangle %d\n",i);
            double x,y,z ;

            Triangle *tri = [[scene triangles] objectAtIndex:i] ;
            x = [[tri vertexAt:0] X];
            y = [[tri vertexAt:0] Y];
            z = [[tri vertexAt:0] Z];
            printf("    A      : %3.6f,%3.6f,%3.6f\n",x,y,z );
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&z, sizeof(double), 1, fp);
            x = [[tri vertexAt:1] X];
            y = [[tri vertexAt:1] Y];
            z = [[tri vertexAt:1] Z];
            printf("    B      : %3.6f,%3.6f,%3.6f\n",x,y,z );
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&z, sizeof(double), 1, fp);
            x = [[tri vertexAt:2] X];
            y = [[tri vertexAt:2] Y];
            z = [[tri vertexAt:2] Z];
            printf("    C      : %3.6f,%3.6f,%3.6f\n",x,y,z );
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&z, sizeof(double), 1, fp);
            
            
            x = [[tri normal] X];
            y = [[tri normal] Y];
            z = [[tri normal] Z];
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&z, sizeof(double), 1, fp);
            printf("    Normal : %3.6f,%3.6f,%3.6f\n",x,y,z);
            
            char matname[MAXMATNAME];
            int matind = i % nmat;
            radarMaterial *material;
            material = [[scene triangleMaterials] objectAtIndex:matind] ;
            strcpy(matname, [[material materialName] UTF8String]);
            fwrite(matname, sizeof(char), MAXMATNAME, fp);
            printf("  Material : %s\n",matname);

        }
        fclose(fp);
        im_close_lib(&status);
        
    }
    return 0;
}

