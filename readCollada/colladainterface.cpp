#include "colladainterface.h"
extern "C"{
#include "matrixMultiplication.h"
}

char array_types[7][15] = {"float_array", "int_array", "bool_array", "Name_array",
    "IDREF_array", "SIDREF_array", "token_array"};

char primitive_types[7][15] = {"lines", "linestrips", "polygons", "polylist",
    "triangles", "trifans", "tristrips"};

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
        TiXmlElement* instance_geometry =
        doc.RootElement()->FirstChildElement("library_visual_scenes")->FirstChildElement("visual_scene")->FirstChildElement("node")->FirstChildElement("instance_geometry");
        while (instance_geometry != NULL) {
            std::string instance_id ;
            instance_id = instance_geometry->Attribute("url") ;
            instance_id = instance_id.erase(0, 1);
            if (data.name == instance_id) {
                std::string materialID ;
                materialID = instance_geometry->FirstChildElement("bind_material")->FirstChildElement("technique_common")->FirstChildElement("instance_material")->Attribute("target");
                materialID = materialID.erase(0, 1);
                std::string materiaName ;
                TiXmlElement *material =  doc.RootElement()->FirstChildElement("library_materials")->FirstChildElement("material");
                while( material != NULL){
                    if(std::string(material->Attribute("id")) == materialID){
                        data.material = std::string(material->Attribute("name")) ;
                    }
                    material = material->NextSiblingElement("material");
                }
            }
            instance_geometry = instance_geometry->NextSiblingElement("instance_geometry");
        }
        
        // Find the world transformation matrix for this geometry
        //
        
        TiXmlElement * node =
        doc.RootElement()->FirstChildElement("library_nodes");
        if(node != NULL){
            printf("Error: Collada file is hierarchical. Please save without \'preserve hierachies\' option");
            exit(1);
        }
        node = doc.RootElement()->FirstChildElement("library_visual_scenes")->FirstChildElement("visual_scene")->FirstChildElement("node");
        double transform[16] = {1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0} ;

        while (node != NULL){
        
            instance_geometry = node->FirstChildElement("instance_geometry");
            double tMat[16];
            
            while (instance_geometry != NULL) {
                std::string instance_id ;
                instance_id = instance_geometry->Attribute("url") ;
                instance_id = instance_id.erase(0, 1);
                if (data.name == instance_id) {
                    // Found the instance_geometry for our geometry
                    // see if it has a transformation matrix
                    //
                    TiXmlElement *matrix = node->FirstChildElement("matrix");
                    if( matrix != NULL){
                        char* matrixText = (char *)matrix->GetText() ;
                        tMat[0] = atof(strtok(matrixText, " "));
                        for(int i=1; i<16; i++){
                            tMat[i] = atof(strtok(NULL, " "));
                        }
                        double outMat[16] ;
                        matmul(transform, tMat, outMat, 4, 4, 4, 4);
                        for(int i=0; i<16; i++)transform[i]=outMat[i];
                    }
                    
                    // Iterate back to top level calculate the transform for each parent node
                    //
                    TiXmlElement *parentMatrix = node->Parent()->FirstChildElement("matrix");
                    while (parentMatrix != NULL){
                        char* matrixText = (char *)parentMatrix->GetText() ;
                        tMat[0] = atof(strtok(matrixText, " "));
                        for(int i=1; i<16; i++){
                            tMat[i] = atof(strtok(NULL, " "));
                        }
                        double outMat[16] ;
                        matmul(transform, tMat, outMat, 4, 4, 4, 4);
                        for(int i=0; i<16; i++)transform[i]=outMat[i];
                        
                        parentMatrix = parentMatrix->Parent()->FirstChildElement("matrix");
                    }
                    // Now have the world transformation matrix stored in transform[]
                    //
                    
                }
                instance_geometry = instance_geometry->NextSiblingElement("instance_geometry");
            }
        
            TiXmlElement *nextnode = node->NextSiblingElement("node");
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
            vertices = mesh->FirstChildElement("vertices");
            input = vertices->FirstChildElement("input");
            
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
            
            while ( primitive != NULL) {
                
                // Determine number of primitives
                //
                primitive->QueryIntAttribute("count", &prim_count);
                num_indices = prim_count * 3;
                data.index_count += num_indices;
                
                // Read the index values
                //
                char* text = (char*)(primitive->FirstChildElement("p")->GetText());
                indices_vec.push_back((unsigned short)atoi(strtok(text, " ")));
                for(int index=1; index<num_indices; index++) {
                    indices_vec.push_back((unsigned short)atoi(strtok(NULL, " ")));
                }
                
                primitive = primitive->NextSiblingElement("triangles");
            }
            
            if(data.index_count > 0){
                data.primitive = GL_TRIANGLES ;
                // Allocate memory for indices
                //
                data.indices = (unsigned short*)malloc(indices_vec.size()* sizeof(unsigned short));
                for (int i=0; i<indices_vec.size(); i++){
                    data.indices[i] = indices_vec[i] ;
                }
            }

            /*
            for(int i=0; i<7; i++) {
                primitive = mesh->FirstChildElement(primitive_types[i]);
                if(primitive != NULL) {
                    
                    // Determine number of primitives
                    //
                    primitive->QueryIntAttribute("count", &prim_count);
                    
                    // Determine primitive type and set count
                    //
                    switch(i) {
                        case 0:
                            data.primitive = GL_LINES;
                            num_indices = prim_count * 2;
                            break;
                        case 1:
                            data.primitive = GL_LINE_STRIP;
                            num_indices = prim_count + 1;
                            break;
                        case 4:
                            data.primitive = GL_TRIANGLES;
                            num_indices = prim_count * 3;
                            break;
                        case 5:
                            data.primitive = GL_TRIANGLE_FAN;
                            num_indices = prim_count + 2;
                            break;
                        case 6:
                            data.primitive = GL_TRIANGLE_STRIP;
                            num_indices = prim_count + 2;
                            break;
                        default: std::cout << "Primitive " << primitive_types[i] <<
                            " not supported" << std::endl;
                    }
                    data.index_count = num_indices;
                    
                    // Allocate memory for indices
                    //
                    data.indices = (unsigned short*)malloc(num_indices * sizeof(unsigned short));
                    
                    // Read the index values
                    //
                    char* text = (char*)(primitive->FirstChildElement("p")->GetText());
                    data.indices[0] = (unsigned short)atoi(strtok(text, " "));
                    for(int index=1; index<num_indices; index++) {
                        data.indices[index] = (unsigned short)atoi(strtok(NULL, " "));
                    }
                }
            }*/
            
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
        free(geom_it->indices);
        
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
    char* text;
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
            text = (char*)(array->GetText());
            
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
                    ((float*)source_data.data)[0] = atof(strtok(text, " "));  
                    for(unsigned int index=1; index<num_vals; index++) {
                        ((float*)source_data.data)[index] = atof(strtok(NULL, " "));   
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
                    ((int*)source_data.data)[0] = atof(strtok(text, " "));  
                    for(unsigned int index=1; index<num_vals; index++) {
                        ((int*)source_data.data)[index] = atof(strtok(NULL, " "));   
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


