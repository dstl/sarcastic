

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "colladainterface.h"
#include <SIlib/SIlib.h>
#include "trianglecpp.h"
#include "materialProperties.h"

std::vector<ColGeom> geom_vec;    // Vector containing COLLADA meshes
int num_objects;                  // Number of meshes in the vector
int triangleCount(ColGeom v);
int vertexCount(ColGeom v);
bool SPVectorsAreEqual(SPVector a, SPVector b);
std::vector<SPVector> uniqueVertices(std::vector<Triangle> triangles);
void writePLYFile(std::string fname, std::vector<Triangle> triangles);
FILE * initialise_ply_file(const char *fname, int nvertices, int ntris);

int main(int argc, char* argv[]) {
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"readCollada", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    char * inFile, * outFile ;
    bool plyop=false;
    char *plyFname;

    inFile = input_string("Collada .dae file to read", "inFile",
                          "The pathname of a Collada file containing the scene",
                          "/Users/Darren/Development/Models/scene.dae");
    
    outFile = input_string("Name of output file", "outFile",
                           "The pathname of the file to create. This will be a binary file containing triangle and material information of the scene",
                           "/tmp/triangles.tri");
    
    plyop = input_yesno("Write out collada file as PLY file", "plyop",
                        "The PLY format is a simple-to-understand representation ofteh transgles in the Collada scene",
                        plyop);
    if(plyop){
        plyFname = input_string("name of PLY file to create", "plyFname",
                               "Full path name of the PLY file to create", plyFname);
    }
    
    
    // Geom_vec is a C++ Vector containing an array of 'meshes' from the collada file
    // We are only interested in the triangles in the file so loop though the
    // Geom_vec array, find all the triangles then write them to a triangle vector
    //
    ColladaInterface::readGeometries(&geom_vec, inFile);
    
    std::vector<Triangle> tri_vec;
    num_objects = (int) geom_vec.size();
    for(int i=0; i<num_objects; i++){
        if(geom_vec[i].primitive == GL_TRIANGLES ){
            float *f      = (float *)geom_vec[i].map["POSITION"].data;
            int stride    = geom_vec[i].map["POSITION"].stride ;
            float scaling = geom_vec[i].scaling;

            for (int j=0; j<geom_vec[i].index_count/3; j++){
                short AAi,BBi,CCi;
                AAi = geom_vec[i].indices[j*3+0] ;
                BBi = geom_vec[i].indices[j*3+1] ;
                CCi = geom_vec[i].indices[j*3+2] ;
                float AAx,AAy,AAz, BBx,BBy,BBz, CCx,CCy,CCz;
                AAx = f[AAi*stride+0] * scaling;
                AAy = f[AAi*stride+1] * scaling;
                AAz = f[AAi*stride+2] * scaling;
                BBx = f[BBi*stride+0] * scaling;
                BBy = f[BBi*stride+1] * scaling;
                BBz = f[BBi*stride+2] * scaling;
                CCx = f[CCi*stride+0] * scaling;
                CCy = f[CCi*stride+1] * scaling;
                CCz = f[CCi*stride+2] * scaling;
                SPVector AA,BB,CC;
                VECT_CREATE(AAx, AAy, AAz, AA);
                VECT_CREATE(BBx, BBy, BBz, BB);
                VECT_CREATE(CCx, CCy, CCz, CC);
                Triangle t = Triangle(AA, BB, CC,geom_vec[i].material);
                tri_vec.push_back(t);
            }
        }
    }
    ColladaInterface::freeGeometries(&geom_vec);
    
    
    writePLYFile(plyFname, tri_vec);
    
    
    
    im_close_lib(&status);

    std::cout << "Done" << std::endl;
    
    return 0;
}

int vertexCount(ColGeom v){
    int nvertices = 0 ;
    if(v.map["POSITION"].type == GL_FLOAT){
        nvertices += (v.map["POSITION"].size / sizeof(float))/3 ;
    }else if(v.map["POSITION"].type == GL_INT){
        nvertices += (v.map["POSITION"].size / sizeof(int))/3 ;
    }else{
        printf("Error : data source is neither float nor int\n");
        exit(1);
    }
    return nvertices ;
}

int triangleCount(ColGeom v){
    int ntris = 0;
    
    if(v.primitive == GL_TRIANGLES){
        ntris += v.index_count/3 ;
    }
    return ntris ;
}

FILE * initialise_ply_file(const char *fname, int nvertices, int ntris){
    
    FILE *fp = fopen(fname, "w");
    
    if (fp == NULL) {
        printf("Error : could not open file for writing\n");
        exit(1);
    }
    
    fprintf(fp,"ply\n");
    fprintf(fp,"format ascii 1.0\n");
    fprintf(fp,"comment testplyfile\n");
    fprintf(fp,"element vertex %d\n",nvertices);
    fprintf(fp,"property float x\n");
    fprintf(fp,"property float y\n");
    fprintf(fp,"property float z\n");
    fprintf(fp,"element face %d\n",ntris);
    fprintf(fp,"property list uchar int vertex_index\n");
    fprintf(fp,"property uchar red\n");
    fprintf(fp,"property uchar green\n");
    fprintf(fp,"property uchar blue\n");
    fprintf(fp,"end_header\n");
    
    return(fp);
}

std::vector<SPVector> uniqueVertices(std::vector<Triangle> triangles){
    int ntris = (int)triangles.size();
    SPVector test;
    bool repeated;
    
    std::vector<SPVector> vertices_orig;
    std::vector<SPVector> vertices_new;
    
    for(int i=0; i<ntris; i++){
        vertices_orig.push_back(triangles[i].AA);
        vertices_orig.push_back(triangles[i].BB);
        vertices_orig.push_back(triangles[i].CC);
    }
    
    for(int i=0; i<vertices_orig.size(); i++){
        test     = vertices_orig[i];
        repeated = false;
        for (int j=0; j<vertices_new.size(); j++){
            if( SPVectorsAreEqual(test, vertices_new[j])){
                repeated = true;
            }
        }
        if (repeated == false) {
            vertices_new.push_back(test);
        }
    }
    return vertices_new ;
}

void writePLYFile(std::string fname, std::vector<Triangle> triangles){
    
    std::vector<SPVector> v;
    
    v = uniqueVertices(triangles);
    
    FILE * fp = initialise_ply_file(fname.c_str(), (int)v.size(), (int)triangles.size());
    
    // rebuild triangles to use unique vertices
    //
    int t[(int)triangles.size()][3];
    for(int i=0; i<triangles.size(); i++){
        SPVector aa,bb,cc;
        aa = triangles[i].AA ;
        bb = triangles[i].BB ;
        cc = triangles[i].CC ;
        
        for(int j=0; j<v.size(); j++){
            if(SPVectorsAreEqual(v[j], aa))t[i][0] = j ;
            if(SPVectorsAreEqual(v[j], bb))t[i][1] = j ;
            if(SPVectorsAreEqual(v[j], cc))t[i][2] = j ;
        }
    }
    
    // Write out the unique vertices
    //
    for (int i=0; i<v.size(); i++){
        fprintf(fp,"%4.4f %4.4f %4.4f\n",v[i].x,v[i].y,v[i].z);
    }
    // now write out triangles
    //
    for(int i=0; i<triangles.size(); i++){
        fprintf(fp, "3 %d %d %d %d %d %d\n",t[i][0], t[i][1], t[i][2],
                materialColours[triangles[i].matId][0],materialColours[triangles[i].matId][1],materialColours[triangles[i].matId][2]);
    }
    fclose(fp);
    return ;
}

bool SPVectorsAreEqual(SPVector a, SPVector b){
    return (a.x==b.x && a.y==b.y && a.z==b.z) ;
}