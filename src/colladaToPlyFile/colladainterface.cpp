/***************************************************************************
 * 
 *           Module :  colladainterface.cpp
 *          Program :  colladaToPlyFile
 *       Created by :  Darren Muff on Wed Jul 11 08:12:36 2018 -0400
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  A basic C++ class for reading collada (.dae) files
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

#include "colladainterface.h"
extern "C"{
#include "matrixMultiplication.h"
}

char array_types[7][15] = {"float_array", "int_array", "bool_array", "Name_array",
    "IDREF_array", "SIDREF_array", "token_array"};

char primitive_types[7][15] = {"lines", "linestrips", "polygons", "polylist",
    "triangles", "trifans", "tristrips"};
std::vector<std::string> stringTok(std::string s) ;
void ColladaInterface::readGeometries(std::vector<ColGeom>* v, const char* filename) {
    
    TiXmlElement *mesh, *vertices, *input, *source, *primitive;
    std::string source_name;
    int prim_count, num_indices= -1;
    float scaling=1.0;
    
    // Create document and load COLLADA file
    //
    TiXmlDocument doc(filename);
    doc.LoadFile();
    
    if (!strcmp(doc.FirstChild()->Value(), "COLLADA")) {
        printf("Error : Not a Collada file\n");
        exit(0);
    }
    
    TiXmlElement* unit = doc.RootElement()->FirstChildElement("asset")->FirstChildElement("unit");
    if(unit != NULL){
        scaling = atof(unit->Attribute("meter"));
    }
    
    TiXmlElement* geometry =
    doc.RootElement()->FirstChildElement("library_geometries")->FirstChildElement("geometry");
    
    double transform[16];
    double tMat[16];
    TiXmlElement *matrix;
    TiXmlElement *nextnode ;
    TiXmlElement *parentMatrix;
    TiXmlElement* instance_geometry ;
    TiXmlElement * node ;
    
    // Iterate through geometry elements
    //
    while(geometry != NULL) {
        
        // Create new geometry
        //
        ColGeom data;
        
        // set the scaling required to make everything in metres
        //
        data.scaling = scaling ;
        
        // Set the geometry name
        //
        data.name = geometry->Attribute("id");
        
        // initialise the number of indices
        //
        data.index_count = 0;
        
        // Find the material for this geometry
        //
        node = doc.RootElement()->FirstChildElement("library_nodes");
        if(node != NULL){
            printf("Error: Collada file is hierarchical. Please save without \'preserve hierachies\' option");
            exit(1);
        }
        node = doc.RootElement()->FirstChildElement("library_visual_scenes")->FirstChildElement("visual_scene")->FirstChildElement("node");
        while (node != NULL){
            instance_geometry = node->FirstChildElement("instance_geometry");

            while (instance_geometry != NULL) {
                std::string instance_id ;
                instance_id = instance_geometry->Attribute("url") ;
                instance_id = instance_id.erase(0, 1);
                if (data.name == instance_id) {
                    // Found the instance_geometry for our geometry
                    // see if it has a material
                    //
                    std::string materialID ;
                    TiXmlElement *instanceMaterial;
                    instanceMaterial = instance_geometry->FirstChildElement("bind_material")->FirstChildElement("technique_common")->FirstChildElement("instance_material");
                    while (instanceMaterial != NULL) {
                    
                        materialID = instanceMaterial->Attribute("target");
                        materialID = materialID.erase(0, 1);
                        std::string materiaName ;
                        TiXmlElement *material =  doc.RootElement()->FirstChildElement("library_materials")->FirstChildElement("material");
                        while (material != NULL) {
                            if(std::string(material->Attribute("id")) == materialID){
                                
                                if(data.materialSide1 == std::string("")){
                                    data.materialSide1 = std::string(material->Attribute("name")) ;
                                } else if (data.materialSide2 == std::string("")) {
                                    data.materialSide2 = std::string(material->Attribute("name")) ;
                                }
                            }
                            material = material->NextSiblingElement("material");
                        }
                        instanceMaterial = instanceMaterial->NextSiblingElement("instance_material");
                    }
                }
                instance_geometry = instance_geometry->NextSiblingElement("instance_geometry");
            }
            
            // Make sure we also search for child nodes
            //
            TiXmlElement *nextnode = node->NextSiblingElement("node");
            if (nextnode == NULL){
                nextnode = node->FirstChildElement("node");
            }
            node = nextnode;
        }
        
        // Find the world transformation matrix for this geometry
        //
        node = doc.RootElement()->FirstChildElement("library_visual_scenes")->FirstChildElement("visual_scene")->FirstChildElement("node");
        transform[0]  = 1.0; transform[1]  = 0.0; transform[2]  = 0.0; transform[3]  = 0.0;
        transform[4]  = 0.0; transform[5]  = 1.0; transform[6]  = 0.0; transform[7]  = 0.0;
        transform[8]  = 0.0; transform[9]  = 0.0; transform[10] = 1.0; transform[11] = 0.0;
        transform[12] = 0.0; transform[13] = 0.0; transform[14] = 0.0; transform[15] = 1.0;

        while (node != NULL){
            
            instance_geometry = node->FirstChildElement("instance_geometry");
            matrix = node->FirstChildElement("matrix");
            
            while (instance_geometry != NULL) {
                std::string instance_id ;
                instance_id = instance_geometry->Attribute("url") ;
                instance_id = instance_id.erase(0, 1);
                if (data.name == instance_id) {
                    if( matrix != NULL){
                        std::vector<std::string> matrixElements = stringTok(std::string(matrix->GetText()));
                        for(int i=0; i<16; i++){
                            tMat[i] = atof(matrixElements[i].c_str());
                        }
                        double outMat[16] ;
                        matmul(transform, tMat, outMat, 4, 4, 4, 4);
                        for(int i=0; i<16; i++)transform[i]=outMat[i];
                    }
                    // Iterate back to top level calculate the transform for each parent node
                    //
                    TiXmlNode *parentNode = node->Parent();
                    while (parentNode != NULL){
                        parentMatrix = parentNode->FirstChildElement("matrix");
                        if (parentMatrix != NULL) {
                            std::vector<std::string> matrixElements = stringTok(std::string(parentMatrix->GetText()));
                            for(int i=0; i<16; i++){
                                tMat[i] = atof(matrixElements[i].c_str());
                            }
                            double outMat[16] ;
                            matmul(transform, tMat, outMat, 4, 4, 4, 4);
                            for(int i=0; i<16; i++)transform[i]=outMat[i];
                        }
                        parentNode   = parentNode->Parent() ;
                    }
                    // Now have the world transformation matrix stored in transform[]
                    //
                }
                instance_geometry = instance_geometry->NextSiblingElement("instance_geometry");
            }
            nextnode = node->NextSiblingElement("node");
            if (nextnode == NULL){
                nextnode = node->FirstChildElement("node");
            }
            node = nextnode;
        }
        
        for(int i=0; i<16; i++)data.transform[i]=transform[i];

        // Iterate through mesh elements
        //
        mesh = geometry->FirstChildElement("mesh");
        while(mesh != NULL) {
            input = NULL;
            vertices = mesh->FirstChildElement("vertices");
            if(vertices != NULL){
                input = vertices->FirstChildElement("input");
            }
            
            // Iterate through input elements
            //
            while(input != NULL) {
                source_name = std::string(input->Attribute("source"));
                source_name = source_name.erase(0, 1);
                source = mesh->FirstChildElement("source");
                
                // Iterate through source elements
                //
                while(source != NULL) {
                    if(std::string(source->Attribute("id")) == source_name) {
                        data.map[std::string(input->Attribute("semantic"))] = readSource(source);
                    }
                    
                    source = source->NextSiblingElement("source");
                }
                
                input = input->NextSiblingElement("input");
            }
            
            // If there are any triangles then read them in
            // first find total number of triangles
            //
            primitive = mesh->FirstChildElement("triangles");
            std::vector<unsigned short>indices_vec ;
            
            if(primitive != NULL) {
                
                // Determine number of primitives
                //
                primitive->QueryIntAttribute("count", &prim_count);
                
                data.primitive = GL_TRIANGLES;
                num_indices = prim_count * 3;
                data.index_count = num_indices;
                
                // Allocate memory for indices
                //
                if(num_indices == 0){
                    printf("Error: Number of indices is zero\n");
                    exit(0);
                }
                data.indices = (unsigned short*)malloc(num_indices * sizeof(unsigned short));
                
                // Read the index values
                //
                std::vector<std::string> text = stringTok(std::string(primitive->FirstChildElement("p")->GetText()));
                
                for(int index=0; index<num_indices; index++) {
                    data.indices[index] = (unsigned short)atoi(text[index].c_str());
                }
            }
            
            mesh = mesh->NextSiblingElement("mesh");
        }
        
        v->push_back(data);
        
        geometry = geometry->NextSiblingElement("geometry");
    }
}

void ColladaInterface::freeGeometries(std::vector<ColGeom>* v) {
    
    std::vector<ColGeom>::iterator geom_it;
    SourceMap::iterator map_it;
    
    for(geom_it = v->begin(); geom_it < v->end(); geom_it++) {
        
        // Deallocate index data
        //
        if(geom_it->index_count > 0){
            free(geom_it->indices);
        }
        
        // Deallocate array data in each map value
        //
        for(map_it = geom_it->map.begin(); map_it != geom_it->map.end(); map_it++) {
            free((*map_it).second.data);
        }
        
        // Erase the current ColGeom from the vector
        //
        v->erase(geom_it);
    }
}

SourceData readSource(TiXmlElement* source) {
    
    SourceData source_data;
    TiXmlElement *array;
    std::vector<std::string> text;
    unsigned int num_vals, stride;
    int check;
    
    for(int i=0; i<7; i++) {
        array = source->FirstChildElement(array_types[i]);
        if(array != NULL) {
            
            // Find number of values
            //
            array->QueryUnsignedAttribute("count", &num_vals);
            source_data.size = num_vals;
            
            // Find stride
            //
            check = source->FirstChildElement("technique_common")->FirstChildElement("accessor")->QueryUnsignedAttribute("stride", &stride);
            if(check != TIXML_NO_ATTRIBUTE)
                source_data.stride = stride;
            else
                source_data.stride = 1;
            
            // Read array values
            //
            text = stringTok(array->GetText());
            
            // Initialize mesh data according to data type
            //
            switch(i) {
                    
                    // Array of floats
                    //
                case 0:
                    source_data.type = GL_FLOAT;
                    source_data.size *= sizeof(float);
                    source_data.data = malloc(num_vals * sizeof(float));
                    
                    // Read the float values
                    //
                    for(unsigned int index=0; index<num_vals; index++) {
                        ((float*)source_data.data)[index] = atof(text[index].c_str());
                    }
                    break;
                    
                    // Array of integers
                    //
                case 1:
                    source_data.type = GL_INT;
                    source_data.size *= sizeof(int);
                    source_data.data = malloc(num_vals * sizeof(int));
                    
                    // Read the int values
                    //
                    for(unsigned int index=0; index<num_vals; index++) {
                        ((int*)source_data.data)[index] = atof(text[index].c_str());
                    }
                    break;
                    
                    // Other
                    //
                default:
                    std::cout << "Collada Reader doesn't support mesh data in this format" << std::endl;
                    break;
            }
        }
    }
    return source_data;
}

std::vector<std::string> stringTok(std::string s){
    std::vector<std::string> results;
    size_t lastOffset = 0;
    std::string delims(" ");
    
    while(true)
    {
        size_t offset = s.find_first_of(delims, lastOffset);
        results.push_back(s.substr(lastOffset, offset - lastOffset));
        if (offset == std::string::npos)
            break;
        else
            lastOffset = offset + 1; // add one to skip the delimiter
    }
    
    return results ;
}
