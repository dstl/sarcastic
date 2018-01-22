/***************************************************************************
 *
 *       Module:    readMaterialFile.hpp
 *      Program:    sarcastic (and others)
 *   Created by:    Darren Muff on 8/09/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *       This programme reads in material properties from a text file.
 *   The material properties file :
 *       Must use the '#' caharacter as the first character of a line to exclude it
 *       Must contain the text : "# Material Properties File" as the first line
 *       Must contain the text " "# Version 1.0" as the second line where ?? is the
 *       version number defined in compatMatFileVersion (in readMaterialFile.hpp)
 *   The file can tain any number of materials but the file must contain one line
 *   of the format "NMATERIALS ??" where ?? is the number of materials in the file.
 *   Each material must be of the format:
 *       int charstring  float   float   float float float float float int int int
 *   corresponding to:
 *       id  Name  corrLen Roughness Rs Rm Specular Diffuse Shinyness R G B
 *   where :
 *       id is the index of this material used to check that the file contains the correct number of materials
 *       Name is a text field describing the material
 *       corrLen is the correlation length of the material in metres (ie over what distance is it flat)
 *       Roughness is the mean surface roughness of the material
 *       Rs is the electrical resistivity of the material in Ohm.Metres
 *       Rm is the magnetic condictivity of the material
 *       Specular is the amount of incident light that reflects in the specular direction
 *       Diffuse is the amount of incident light that reflects diffusely.
 *       Shinyness is the alpha term in the phong shading equation and defines the solid angle that
 *       the specular scattering direction covers. (Higher values are smaller, more specular).
 *       R G B are used to define the colour of the material in any output files.
 *
 *   CLASSIFICATION        :  OFFICIAL
 *   Date of CLASSN        :  8/09/2017
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
 *
 ***************************************************************************/

#ifndef readMaterialFile_hpp
#define readMaterialFile_hpp
#define compatMatFileVersion "1.0"

#include <iostream>
#include "materialProperties.h"
#include <fstream>
#include <sstream>
#include <regex>

/// Function to read in material properties from a file
/// fileName is input
/// nMaterials is return - the number of materials in the file
/// matProperties is returned. A sctruture containing material information
/// matColours is returned. Colours used for teh different materials
/// Memeory is allocated for matProperties and matColours and it is the responsibility of teh calling function
/// to free up this memory
///
void readMaterialFile(std::string fileName, bool printOutput) ;
void initialiseMaterials(std::string filename, bool printOutput);

#endif /* readMaterialFile_hpp */
