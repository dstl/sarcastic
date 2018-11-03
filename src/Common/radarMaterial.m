/***************************************************************************
 * 
 *           Module :  radarMaterial.m
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 19/05/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description : Objective C class for radar materials
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

#import "radarMaterial.h"

@implementation radarMaterial
@synthesize materialName ;
@synthesize materialID ;
@synthesize specular ;
@synthesize diffuse ;
@synthesize shinyness ;
@synthesize ambient ;
@synthesize roughness ;
@synthesize correlationPatchSize ;
@synthesize materialTypeArray ;
- (id) initWithMaterialName: (NSString *) name {
    self = [super init];
    if (self) {
        materialTypeArray = [NSArray arrayWithObjects:
            @"asphalt",
            @"brick",
            @"concrete",
            @"metal",
            @"material",
            @"roofing",
            @"vegetation",
            @"water",
            @"wood", nil];

        materialID   = [self materialTypeStringToEnum:name] ;
        materialName = [self materialTypeEnumToString:materialID];
        ambient = 0;
        
        switch (materialID) {
            case ASPHALT :
                specular  = 0.8f;
                diffuse   = 0.2f;
                shinyness = 30.0f;
                break;
            case BRICK :
                specular  = 0.7f;
                diffuse   = 0.3f;
                shinyness = 20.0f;
                break;
            case CONCRETE :
                specular  = 0.3f;
                diffuse   = 0.7f;
                shinyness = 10.0f;
                break;
            case METAL :
                specular  = 0.1f;
                diffuse   = 0.0f;
                shinyness = 30.0f;
                break;
            case ROOFING :
                specular  = 0.7f;
                diffuse   = 0.3f;
                shinyness = 40.0f;
                break;
            case VEGETATION :
                specular  = 0.2f;
                diffuse   = 0.8f;
                shinyness = 5.0f;
                break;
            case WATER :
                specular  = 1.0f;
                diffuse   = 0.0f;
                shinyness = 50.0;
                break;
            case WOOD :
                specular  = 0.6f;
                diffuse   = 0.4f;
                shinyness = 10.0;
                break;
            case MATERIAL :
            default:        // default MATERIAL
                specular  = 0.9f;
                diffuse   = 0.0f;
                shinyness = 30.0;
                break;
        }
    }
    return self;
}

+ (id) radarMaterialWithMaterialName:(NSString *)name{
    return ([[radarMaterial alloc] initWithMaterialName:name]);
}

- (NSString *) materialTypeEnumToString: (materialType) enumVal {
    return materialTypeArray[enumVal];
}

- (materialType) materialTypeStringToEnum: (NSString *) strVal {
    NSString *material;
    NSRange range ;
    for (material in materialTypeArray){
        range = [strVal rangeOfString:material options:NSCaseInsensitiveSearch];
        if (range.location != NSNotFound) {
            return ((int)[materialTypeArray indexOfObject:material]);
        }
    }
    // Got here but no solution so use defaut type
    //
    for (material in materialTypeArray){
        range = [material rangeOfString:@"material" options:NSCaseInsensitiveSearch];
        
        if (range.location != NSNotFound) {
            return ((int)[materialTypeArray indexOfObject:material]);
        }
    }
    return 0;

}
@end
