/***************************************************************************
 *  COLScene.m
 *  Sadilac
 *
 *  Created by Muff Darren on 19/05/2012.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   18/02/2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
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
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 ***************************************************************************/

#import "COLScene.h"

@implementation COLScene 
@synthesize rootnode=_rootnode ;
@synthesize triangles=_triangles;
@synthesize triangleMaterials=_triangleMaterials;

- (id) initWithURL:(NSURL *)url {
    self = [super init] ;
    if (self) {
        NSError * __autoreleasing error;
        SCNScene *scene = [SCNScene sceneWithURL:url options:nil error:&error];
        if (!scene){
            NSLog(@"DAE file %s not found!",[[url absoluteString] UTF8String]);
            exit(1);
        }
        
        // Initialiese the triangle array
        //
        _triangles = [[NSMutableOrderedSet alloc] init];
        
        // Assign rootnode so we can get it later
        //
        _rootnode     = scene.rootNode ;
        
        // Work out the sclaing used within the collada file
        //
        [self readScalingFromXML:url];
        
        // set up an array of materials we know about
        //
        _triangleMaterials = [[NSMutableOrderedSet alloc] init];
        
        // Build the triangle array from the collada scene primitives
        //
        [self buildTriangles];
        
    }
    return self ;
}

+ (id) sceneWithURL: (NSURL *) url { return [[COLScene alloc] initWithURL: url]; }


- (void)readScalingFromXML:(NSURL *)url {
    
    NSData *data = [NSData dataWithContentsOfURL:url];
    NSXMLParser* parser = [[NSXMLParser alloc] initWithData: data];
    
    [parser setDelegate:self];
    [parser parse];
}

- (void)parser:(NSXMLParser *)parser didStartElement:(NSString *)elementName namespaceURI:(NSString *)namespaceURI qualifiedName:(NSString *)qualifiedName attributes:(NSDictionary *)attributeDict{
    
    if ([elementName isEqualToString:@"unit"]) {
        
        float scaling = [[attributeDict valueForKey:@"meter"] floatValue];
//        NSString* name = [attributeDict valueForKey:@"name"];
//        NSLog(@"units are: %@, scaling: %f", name, scaling);
        _sceneScaling = scaling;
    }
}

- (void) buildTrianglesForNode: (SCNNode *) node {
    SCNVector3  *vertices = NULL;
    SCNVector3  *normals  = NULL;
    SCNGeometry *geometry = [node geometry];
    
    // Set up the scaling and transformation matrices for this node
    //
    CATransform3D scaleMat = CATransform3DMakeScale(_sceneScaling, _sceneScaling, _sceneScaling);
    CATransform3D transMat = CATransform3DConcat([node worldTransform], scaleMat);
    
    // Get the Vertex source data for this geometry description
    //
    NSArray *vertexSources = [geometry geometrySourcesForSemantic:SCNGeometrySourceSemanticVertex];
    
    // Get the source data for the normals
    //
    NSArray *normalSources = [geometry geometrySourcesForSemantic:SCNGeometrySourceSemanticNormal];
    
    // Unpack the vertex data into an array of vertices. (each vertex is of type SCNVector3)
    //
    for (int i=0; i<[vertexSources count]; i++){
        SCNGeometrySource *gs                        = [vertexSources objectAtIndex:i];
                           vertices                  = (SCNVector3 *)malloc(sizeof(SCNVector3)*gs.vectorCount);
        NSInteger          stride                    = gs.dataStride; // in bytes
        NSInteger          offset                    = gs.dataOffset; // in bytes
        NSInteger          componentsPerVector       = gs.componentsPerVector;
        NSInteger          bytesPerVector            = componentsPerVector * gs.bytesPerComponent;
        
        for (int j=0; j<gs.vectorCount; j++) {
            float vectorData[componentsPerVector];
            NSRange byteRange = NSMakeRange(j*stride + offset, bytesPerVector);
            [gs.data getBytes:&vectorData range:byteRange];
            vertices[j]       = SCNVector3Make(vectorData[0], vectorData[1], vectorData[2]);
        }
    }
    
    // Unpack the normal data into an array of normals. (each normal is of type SCNVector3)
    //
    for (int i=0; i<[normalSources count]; i++){
        SCNGeometrySource *gs                        = [normalSources objectAtIndex:i];
        normals                                      = (SCNVector3 *)malloc(sizeof(SCNVector3)*gs.vectorCount);
        NSInteger          stride                    = gs.dataStride; // in bytes
        NSInteger          offset                    = gs.dataOffset; // in bytes
        NSInteger          componentsPerVector       = gs.componentsPerVector;
        NSInteger          bytesPerVector            = componentsPerVector * gs.bytesPerComponent;
        
        for (int j=0; j<gs.vectorCount; j++) {
            float vectorData[componentsPerVector];
            NSRange byteRange = NSMakeRange(j*stride + offset, bytesPerVector);
            [gs.data getBytes:&vectorData range:byteRange];
            normals[j]        = SCNVector3Make(vectorData[0], vectorData[1], vectorData[2]);
        }
    }
    
    // If this geomtry has elements associated with it then find out what type they are and unpack them
    // Each primitive has a set of indices that are used to find the vertex points (from vertices[]) for
    // the primitive
    //
    if ([geometry geometryElementCount] >0) {
        for (int i=0; i<[geometry geometryElementCount]; i++){
            SCNGeometryElement *element   = [geometry geometryElementAtIndex:i];
            NSArray   *materials          = [geometry materials];
            NSInteger nPrimitives         = [element primitiveCount];
            NSInteger bytesPerIndex       = [element bytesPerIndex];
            switch (element.primitiveType) {
                case SCNGeometryPrimitiveTypeTriangles:     // Triangles
                    for (int j=0; j<nPrimitives; j++) {
                        Triangle *tri ;
                        if (bytesPerIndex == 1) {
                            
                            unsigned char *indices;
                            indices = malloc(3*bytesPerIndex);
                            NSRange byteRange = NSMakeRange(j*3*bytesPerIndex, 3*bytesPerIndex);
                            [element.data getBytes:indices range:byteRange];
                            
                            // Now indices are unpacked we need to apply scaling and translation. This is because
                            // Each node in the collada tree can contain primitives which are later translated
                            // relative to the parent node
                            //
                            CATransform3D ans;
                            ans = CATransform3DTranslate(transMat, vertices[indices[0]].x, vertices[indices[0]].y, vertices[indices[0]].z);
                            MVector *a = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            ans = CATransform3DTranslate(transMat, vertices[indices[1]].x, vertices[indices[1]].y, vertices[indices[1]].z);
                            MVector *b = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            ans = CATransform3DTranslate(transMat, vertices[indices[2]].x, vertices[indices[2]].y, vertices[indices[2]].z);
                            MVector *c = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            
                            ans = CATransform3DTranslate(transMat, normals[indices[0]].x, normals[indices[0]].y, normals[indices[0]].z);
                            MVector *norma = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            ans = CATransform3DTranslate(transMat, normals[indices[1]].x, normals[indices[1]].y, normals[indices[1]].z);
                            MVector *normb = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            ans = CATransform3DTranslate(transMat, normals[indices[2]].x, normals[indices[2]].y, normals[indices[2]].z);
                            MVector *normc = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            
                            MVector *normal =[MVector MVectorWithValuesX:([norma X]+[normb X]+[normc X])/3.0
                                                                       Y:([norma Y]+[normb Y]+[normc Y])/3.0
                                                                       Z:([norma Z]+[normb Z]+[normc Z])/3.0] ;
                            [normal normalise];

                            tri = [Triangle TriangleWithVerticesAa:a Bb:b Cc:c andNormal: normal ];
                            
                            free(indices);
                            
                        }else if (bytesPerIndex == 2){
                            
                            unsigned short *indices;
                            indices = malloc(3*bytesPerIndex);
                            NSRange byteRange = NSMakeRange(j*3*bytesPerIndex, 3*bytesPerIndex);
                            [element.data getBytes:indices range:byteRange];
                            
                            // Now indices are unpacked we need to apply scaling and translation. This is because
                            // Each node in the collada tree can contain primitives which are latere translated
                            // relative to the parent node
                            //
                            CATransform3D ans;
                            ans = CATransform3DTranslate(transMat, vertices[indices[0]].x, vertices[indices[0]].y, vertices[indices[0]].z);
                            MVector *a = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            ans = CATransform3DTranslate(transMat, vertices[indices[1]].x, vertices[indices[1]].y, vertices[indices[1]].z);
                            MVector *b = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            ans = CATransform3DTranslate(transMat, vertices[indices[2]].x, vertices[indices[2]].y, vertices[indices[2]].z);
                            MVector *c = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43];
                            
                            ans = CATransform3DTranslate(transMat, normals[indices[0]].x, normals[indices[0]].y, normals[indices[0]].z);
                            MVector *norma = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            ans = CATransform3DTranslate(transMat, normals[indices[1]].x, normals[indices[1]].y, normals[indices[1]].z);
                            MVector *normb = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            ans = CATransform3DTranslate(transMat, normals[indices[2]].x, normals[indices[2]].y, normals[indices[2]].z);
                            MVector *normc = [MVector MVectorWithValuesX:ans.m41 Y:ans.m42 Z:ans.m43] ;
                            
                            MVector *normal =[MVector MVectorWithValuesX:([norma X]+[normb X]+[normc X])/3.0
                                                                       Y:([norma Y]+[normb Y]+[normc Y])/3.0
                                                                       Z:([norma Z]+[normb Z]+[normc Z])/3.0] ;
                            [normal normalise];

                            tri = [Triangle TriangleWithVerticesAa:a Bb:b Cc:c andNormal: normal ];
                            
                            free(indices);

                        }else{
                            printf("Cannot handle bytesPerIndex of %ld\n",(long)bytesPerIndex);
                            exit(0);
                        }
                        
                        // From p44 of Scene Kit Framework Reference, 2014-02-11
                        //  'Each geometry element can be rendered using a different material. The index of the material used for a geometry
                        //   element is equal to the index of that element modulo the number of materials.'
                        //
                        radarMaterial *mat ;
                        
                        mat = [radarMaterial radarMaterialWithMaterialName:[[materials objectAtIndex:0] name]] ;

                        if ([materials count] > 1) {
                            for (int i=0; i< [materials count]; i++){
                                if ( ! [[[materials objectAtIndex:i] name] caseInsensitiveCompare:@"material"] ) {
                                    mat = [radarMaterial radarMaterialWithMaterialName:[[materials objectAtIndex:0] name]] ;
                                }
                            }
                        }
                        [tri setMaterialName: [[materials objectAtIndex:0] name]];
                        [_triangles         addObject:tri];
                        
                        // Check to see if material is already in the list
                        //
                        radarMaterial * m;
                        for ( m in _triangleMaterials ) {
                            if ( [mat materialID] == [m materialID] ) {
                                mat = m ;
                            }
                        }
                        [_triangleMaterials addObject:mat];
                    }
                    break;
                    
                case SCNGeometryPrimitiveTypeTriangleStrip :
                    printf("Ignoring Triangle Strip primitive\n");
                    break;
                case SCNGeometryPrimitiveTypeLine:
                    printf("Ignoring Line primitive\n");
                    break ;
                case SCNGeometryPrimitiveTypePoint :
                    printf("Ignoring point primitive\n");
                    break ;
                default:
                    printf("Primitives other than triangles not supported!\n");
                    break;
            }
        }
    }
    
    free(vertices);
    free(normals);
    
}

- (void) buildTriangles {
    [self traverseTree:_rootnode withSelector:@selector(buildTrianglesForNode:)];
}

- (void) printNode:(SCNNode *) node {
    
    SCNVector3 boxMin, boxMax;
    
    NSString *name = [node name];
    printf("-=-=-=-=--=-=--=-=-=-=-=-=-=-=-=-=-=-=-\n");
    if (name) {
        printf("Node name is %s\n",[name UTF8String]);
    }
    
    printf("Node parent is %s\n",[[[node parentNode] name] UTF8String]);
    
    [node getBoundingBoxMin:&boxMin max:&boxMax];
    float s = _sceneScaling ;
    printf("Box min is %f,%f,%f\nBox Max is %f,%f,%f\n",s*boxMin.x,s*boxMin.y,s*boxMin.z,s*boxMax.x,s*boxMax.y,s*boxMax.z);
    
    SCNCamera *camera = [node camera];
    if(camera)
        printf("Node has a camera\n");
    
    NSArray *children = [node childNodes];
    printf("Node has %ld children\n", [children count]);
    
    // Transformation data
    CATransform3D pivot = [node pivot];
    printf("Pivot matrix:\n");
    printf("  |  %f  %f  %f  %f |\n",s*pivot.m11,s*pivot.m12,s*pivot.m13,s*pivot.m14);
    printf("  |  %f  %f  %f  %f |\n",s*pivot.m21,s*pivot.m22,s*pivot.m23,s*pivot.m24);
    printf("  |  %f  %f  %f  %f |\n",s*pivot.m31,s*pivot.m32,s*pivot.m33,s*pivot.m34);
    printf("  |  %f  %f  %f  %f |\n",s*pivot.m41,s*pivot.m42,s*pivot.m43,s*pivot.m44);
    
    SCNVector3 position = [node position];
    printf("Node position:\n");
    printf("  %f,%f,%f\n",s*position.x,s*position.y,s*position.z);
    
    SCNVector4 rotation = [node rotation];
    printf("Node rotation:\n");
    printf("  axis : %f,%f,%f angle : %f\n",rotation.x,rotation.y,rotation.z,rotation.w);
    
    SCNVector3 scale = [node scale];
    printf("Node Scale:\n");
    printf("  %f,%f,%f\n",scale.x,scale.y,scale.z);
    
    CATransform3D transform = [node transform];
    printf("transform matrix:\n");
    printf("  |  %f  %f  %f  %f |\n",s*transform.m11,s*transform.m12,s*transform.m13,s*transform.m14);
    printf("  |  %f  %f  %f  %f |\n",s*transform.m21,s*transform.m22,s*transform.m23,s*transform.m24);
    printf("  |  %f  %f  %f  %f |\n",s*transform.m31,s*transform.m32,s*transform.m33,s*transform.m34);
    printf("  |  %f  %f  %f  %f |\n",s*transform.m41,s*transform.m42,s*transform.m43,s*transform.m44);
    
    CATransform3D worldTransform = [node worldTransform];
    printf("World Transform matrix:\n");
    printf("  |  %f  %f  %f  %f |\n",worldTransform.m11,worldTransform.m12,worldTransform.m13,worldTransform.m14);
    printf("  |  %f  %f  %f  %f |\n",worldTransform.m21,worldTransform.m22,worldTransform.m23,worldTransform.m24);
    printf("  |  %f  %f  %f  %f |\n",worldTransform.m31,worldTransform.m32,worldTransform.m33,worldTransform.m34);
    printf("  |  %f  %f  %f  %f |\n",worldTransform.m41,worldTransform.m42,worldTransform.m43,worldTransform.m44);
    
    CATransform3D scaleMat = CATransform3DMakeScale(s, s, s);
    CATransform3D transMat = CATransform3DConcat([node worldTransform], scaleMat);
    
    // Geomtery details - first the geometry source data then the geometry elements
    //
    SCNGeometry *geometry = [node geometry];
    if (geometry) {
        printf("Node Geometry details:\n");
        printf("  name is %s\n",[[geometry name] UTF8String]);
        NSArray *vertexSources = [geometry geometrySourcesForSemantic:SCNGeometrySourceSemanticVertex];
        [geometry geometrySourcesForSemantic:SCNGeometrySourceSemanticTexcoord];
        printf("  Geometry Sources:\n");
        printf("    %ld Vertex Sources :\n",[vertexSources count]);
        for (int i=0; i<[vertexSources count]; i++){
            printf("    Vertex source %d:\n",i);
            SCNGeometrySource *gs = [vertexSources objectAtIndex:i];
            printf("      Vectors               : %ld\n",[gs vectorCount]);
            printf("      Components per vector : %ld\n",[gs componentsPerVector]);
            printf("      Bytes per component   : %ld\n",[gs bytesPerComponent]);
            printf("      Data offset in bytes  : %ld\n",[gs dataOffset]);
            printf("      Data stride in bytes  : %ld\n",[gs dataStride]);
            SCNVector3 vertices[gs.vectorCount]; // A new array for vertices
            NSInteger stride = gs.dataStride; // in bytes
            NSInteger offset = gs.dataOffset; // in bytes
            NSInteger componentsPerVector = gs.componentsPerVector;
            NSInteger bytesPerVector = componentsPerVector * gs.bytesPerComponent;
            
            for (int j=0; j<gs.vectorCount; j++) {
                float vectorData[componentsPerVector];
                NSRange byteRange = NSMakeRange(j*stride + offset, bytesPerVector);
                [gs.data getBytes:&vectorData range:byteRange];
                printf("        %f,%f,%f\n",s*vectorData[0],s*vectorData[1],s*vectorData[2]);
                vertices[j] = SCNVector3Make(vectorData[0], vectorData[1], vectorData[2]);
            }
            
            printf("  geometry has %ld elements\n",[geometry geometryElementCount]);
            if ([geometry geometryElementCount] >0) {
                
                NSArray *materials = [geometry materials];

                for (int i=0; i<[geometry geometryElementCount]; i++){
                    SCNGeometryElement *element = [geometry geometryElementAtIndex:i];
                    NSInteger nPrimitives       = [element  primitiveCount];
                    NSInteger bytesPerIndex     = [element  bytesPerIndex];
                    printf("    %ld primitives of type %ld in this element with %ld bytes per index:\n",(long)nPrimitives,[element primitiveType],bytesPerIndex);
                    switch (element.primitiveType) {
                        case 0:     // Triangles
                            for (int tmpi=0; tmpi<nPrimitives; tmpi++) {
                                unsigned char *indices;
                                indices = malloc(3*bytesPerIndex);
                                NSRange byteRange = NSMakeRange(i*3*bytesPerIndex, 3*bytesPerIndex);
                                [element.data getBytes:indices range:byteRange];
                                printf("      %2d  %3d,%3d,%3d\n",tmpi,indices[0*bytesPerIndex],indices[1*bytesPerIndex],indices[2*bytesPerIndex]);
                                free(indices);
                            }
                            printf("    Indices make the following triangles:\n");
                            for (int tmpi=0; tmpi<nPrimitives; tmpi++) {
                                unsigned char *indices;
                                indices = malloc(3*bytesPerIndex);
                                NSRange byteRange = NSMakeRange(tmpi*3*bytesPerIndex, 3*bytesPerIndex);
                                [element.data getBytes:indices range:byteRange];
                                printf("    %d:\n",tmpi);
                                printf("      %f,%f,%f\n",vertices[indices[0*bytesPerIndex]].x,vertices[indices[0*bytesPerIndex]].y,vertices[indices[0*bytesPerIndex]].z);
                                printf("      %f,%f,%f\n",vertices[indices[1*bytesPerIndex]].x,vertices[indices[1*bytesPerIndex]].y,vertices[indices[1*bytesPerIndex]].z);
                                printf("      %f,%f,%f\n",vertices[indices[2*bytesPerIndex]].x,vertices[indices[2*bytesPerIndex]].y,vertices[indices[2]*bytesPerIndex].z);
                                free(indices);
                            }
                            printf("    Triangles after applying scaling and translation matrix:\n");
                            for (int tmpi=0; tmpi<nPrimitives; tmpi++) {
                                unsigned char *indices;
                                indices = malloc(3*bytesPerIndex);
                                NSRange byteRange = NSMakeRange(tmpi*3*bytesPerIndex, 3*bytesPerIndex);
                                [element.data getBytes:indices range:byteRange];
                                printf("    %d:\n",tmpi);
                                for (int j=0; j<3; j++) {
                                    CATransform3D ans = CATransform3DTranslate(transMat, vertices[indices[j*bytesPerIndex]].x, vertices[indices[j*bytesPerIndex]].y, vertices[indices[j*bytesPerIndex]].z);
                                    printf("      %f,%f,%f\n",ans.m41,ans.m42,ans.m43);
                                }
                                free(indices);
                                int matind = tmpi % [materials count];
                                printf("   Triangle material : %s (ind: %d)\n", [[[materials objectAtIndex:matind] name] UTF8String], matind);
                            }
                            break;
                            
                        default:
                            printf("Primitives other than triangles not supported!\n");
                            break;
                    }
                    
                    
                }
                
                if(materials){
                    printf("  geometry has %ld materials\n",[materials count]);
                    SCNMaterial *material;
                    for(int j=0; j<[materials count]; j++){
                        printf("    Material %d\n",j);
                        material = [materials objectAtIndex:j];
                        printf("      Name: %s\n",[[material name] UTF8String]);
                        NSLog(@"      diffuse   : %@",[[material diffuse] contents]);
                        NSLog(@"      specular  : %@",[[material specular] contents]);
                        NSLog(@"      ambient   : %@", [[material ambient] contents]);
                        NSLog(@"      shinyness : %f",[material shininess]);
                      
                        
                    }
                    
                }
            }
        }
    }
    
    SCNLight *light = [node light];
    if (light) {
        printf("Node light details:\n");
    }
    
    return;
}

- (void) printTree {
    [self traverseTree:_rootnode withSelector:@selector(printNode:)] ;
}

// Method that traverses the tree and for each node in the tree it performs
// the method passed in with selector. The method receives a single node object
//
- (void) traverseTree:(SCNNode *) node withSelector: (SEL) selector {
    
    NSArray *children = [node childNodes];
    SCNNode *child;
    if([children count] != 0){
        for (child in children){
            [self traverseTree:child withSelector:selector];
        }
    }
    
    // The original line here was
    //    [self performSelector:selector withObject:node];
    // But this caused a warning as ARC needs to know what to do with the method being called
    // For a full explanations see
    // http://stackoverflow.com/questions/7017281/performselector-may-cause-a-leak-because-its-selector-is-unknown
    //
    IMP imp = [self methodForSelector:selector];
    void (*func)(id, SEL, SCNNode*) = (void *)imp;
    func(self, selector, node);
    
    return ;
}

// Method that parses the name of the texture and returns it as an integer Id
//
//- (int) textureIdForName:(NSString *)textureName{
//    
//    NSString *texture;
//    NSRange range ;
//    for (texture in _triangleMaterials){
//        range = [textureName rangeOfString:texture options:NSCaseInsensitiveSearch];
//        
//        if (range.location != NSNotFound) {
//            return ((int)[_triangleMaterials indexOfObject:texture]);
//        }
//    }
//    return 0;
//}

// Method to print out the triangles in a convenent format to load into grapher
//
- (void) grapherDump {
    int material ;
    Triangle *tri;
    for (material=0; material<[[self triangleMaterials] count]; material++){
        printf("%s\n",[[[self triangleMaterials] objectAtIndex:material] UTF8String]);
        printf("-------------------------------\n");
        for(int i=0; i<[[self triangles] count]; i++){
            tri = [[self triangles] objectAtIndex:i];
            if (!strcmp([[tri materialName] UTF8String], [[[self triangleMaterials] objectAtIndex:material] UTF8String]) ) {
                printf("    %f,%f,%f\n",[[tri vertexAt:0] X],[[tri vertexAt:0] Y],[[tri vertexAt:0] Z]);
                printf("    %f,%f,%f\n",[[tri vertexAt:1] X],[[tri vertexAt:1] Y],[[tri vertexAt:1] Z]);
                printf("    %f,%f,%f\n",[[tri vertexAt:2] X],[[tri vertexAt:2] Y],[[tri vertexAt:2] Z]);
                printf("    %f,%f,%f\n",[[tri vertexAt:0] X],[[tri vertexAt:0] Y],[[tri vertexAt:0] Z]);
            }
        }
    }
}

@end



