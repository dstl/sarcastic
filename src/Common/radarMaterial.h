/***************************************************************************
 * 
 *           Module :  radarMaterial.h
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

#import <Foundation/Foundation.h>

typedef enum {
    ASPHALT,
    BRICK,
    CONCRETE,
    METAL,
    MATERIAL,
    ROOFING,
    VEGETATION,
    WATER,
    WOOD
} materialType ;

@interface radarMaterial : NSObject
@property NSString     * materialName;
@property materialType   materialID;
@property NSArray      * materialTypeArray ;
@property float          specular;
@property float          diffuse;
@property float          shinyness;
@property float          ambient;
@property float          roughness ;
@property float          correlationPatchSize;
- (id) initWithMaterialName: (NSString *) name ;
+ (id) radarMaterialWithMaterialName: (NSString *) name ;
- (NSString *) materialTypeEnumToString: (materialType) enumVal ;
- (materialType) materialTypeStringToEnum: (NSString *) strVal ;

@end
