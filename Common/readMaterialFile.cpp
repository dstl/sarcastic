/***************************************************************************
 *
 *       Module:    readMaterialFile.cpp
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

#include "readMaterialFile.hpp"
using namespace std;

int         globalNumMats ;
scatProps   *globalMatProps;
int         *globalMatColours;

void readMaterialFile(std::string fileName, bool printOutput)
{
    // Validate file
    //
    ifstream file (fileName) ;
    string line ;
    getline(file, line);
    if (line != string("# Material Properties File") ){
        cerr << "Error : input material file does not contain material properties in correct format" << endl;
        exit(1);
    }
    getline(file, line);
    if (line != string("# Version ")+compatMatFileVersion){
        cerr << "Error : Material file version number mismatch "<<endl;
        cerr << "        Material file    : "<<line<<endl;
        cerr << "        Required version : "<<compatMatFileVersion<<endl;
        
    }
    int nl = 0;
    bool nmatsfound = false ;
    while (getline(file, line)) {
        if (line[0] != '#') {
            nl++;
            istringstream iss(line) ;
            string subs;
            iss >> subs ;
            if (subs == string("NMATERIALS")) {
                if (nmatsfound) {
                    cerr << "Number of materials specified twice in materials file"<<endl;
                    exit(1);
                }else{
                    iss >> subs ;
                    globalNumMats = stoi(subs);
                    nmatsfound = true ;
                }
            }
        }
    }
    
    unsigned char *matverify = new unsigned char [globalNumMats] ;
    for(int i=0; i<globalNumMats; i++)matverify[i]=0;
    
    file.clear();
    file.seekg(0) ;
    while (getline(file, line)) {
        if (line[0] != '#') {
            istringstream iss(line) ;
            if( regex_search(line,regex("^\\s*\\d+") )){
                string subs;
                iss >> subs ;
                int ind = stoi(subs) ;
                if(ind >= globalNumMats){
                    cerr << "Inconsistent NMATERIALS number : "<<globalNumMats<<endl;
                    cerr << "Found more than "<<ind<<" materials."<<endl;
                    exit(1);
                }
                matverify[ind] = 1;
            }
        }
    }
    delete [] matverify ;
    
    for(int i=0; i<globalNumMats; i++){
        if ((int)matverify[i] != 1) {
            cerr << "Error : Material id "<<i<<" is not in material file"<<endl;
            exit(1);
        }
    }
        
    globalMatProps   = new scatProps [globalNumMats] ;
    globalMatColours = new int [globalNumMats*3] ;
        
    file.clear() ;
    file.seekg(0) ;
    while (getline(file, line)) {
        if (line[0] != '#') {
            
            istringstream iss(line) ;
            if( regex_search(line,regex("^\\s*\\d+") )){
                
                istringstream testiss(line) ; int testcnt = 0;
                while(testiss){ string substr ; testiss >> substr; testcnt++; }
                if (testcnt != 13){
                    cerr << "Line in material file does not have correct number of fields ("<<testcnt<<")"<<endl;
                    cerr << " Line is : "<<line<<endl;
                    exit(1);
                }
                
                string subs;
                iss >> subs ;
                int ind = stoi(subs) ;
                iss >> subs ;
                strcpy(globalMatProps[ind].matname, subs.c_str());
                iss >> subs ;
                globalMatProps[ind].corlen = stof(subs);
                iss >> subs ;
                globalMatProps[ind].roughness = stof(subs) ;
                iss >> subs ;
                globalMatProps[ind].Rs = stof(subs);
                iss >> subs ;
                globalMatProps[ind].Rm = stof(subs) ;
                iss >> subs ;
                globalMatProps[ind].specular = stof(subs) ;
                iss >> subs ;
                globalMatProps[ind].diffuse = stof(subs) ;
                iss >> subs ;
                globalMatProps[ind].shinyness = stof(subs);
                iss >> subs ;
                globalMatColours[ind*3+0] = stoi(subs) ;
                iss >> subs ;
                globalMatColours[ind*3+1] = stoi(subs) ;
                iss >> subs ;
                globalMatColours[ind*3+2] = stoi(subs) ;
            }
        }
    }
    
    if (printOutput) {
        printf("Material properties used in this file : \n");
        printf("    id         Name CorrLen Roughness    Rs      Rm   Specular Diffuse Shinyness  R   G   B\n");
        for(int i=0; i<globalNumMats; ++i){
            printf("\t%2d ",i);
            printf("%12s ",globalMatProps[i].matname);
            printf("%7.3f ",globalMatProps[i].corlen);
            printf("%9.3f ",globalMatProps[i].roughness);
            printf("%3.1e ",globalMatProps[i].Rs);
            printf("%3.1e ",globalMatProps[i].Rm);
            printf("%8.1f ",globalMatProps[i].specular);
            printf("%7.1f ",globalMatProps[i].diffuse);
            printf("%9.1f ",globalMatProps[i].shinyness);
            printf("%03d ",globalMatColours[i*3+0]);
            printf("%03d ",globalMatColours[i*3+2]);
            printf("%03d ",globalMatColours[i*3+2]);
            printf("\n");
        }
    }
    
    
    return ;
}

void initialiseMaterials(std::string materialsFile, bool printMaterials){
    
    if(string(materialsFile) == string("none") || string(materialsFile) == string("None") || string(materialsFile) == string("")){
        globalNumMats    = NMATERIALS ;
        globalMatProps   = new scatProps [globalNumMats] ;
        globalMatColours = new int [globalNumMats * 3] ;
        for(int i=0; i<globalNumMats; ++i){
            globalMatProps[i]  = materialProperties[i] ;
            globalMatColours[i*3+0] = materialColours[i][0] ;
            globalMatColours[i*3+1] = materialColours[i][1] ;
            globalMatColours[i*3+2] = materialColours[i][2] ;
        }
        
    }else{
        readMaterialFile(string(materialsFile), printMaterials) ;
    }
    
    return ;
}
