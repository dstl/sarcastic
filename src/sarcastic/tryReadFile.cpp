/***************************************************************************
 * 
 *           Module :  tryReadFile.cpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 15/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to repeatedly try and read a file in until the user gets
 *      it right
 *
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

#include "tryReadFile.hpp"

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def)
/// Function to get user input checking that the input is a valid file for reading
///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  key     : a unique bit of text (such as "AzPixelSize") which the input system will then use to store parameter values
///  help    : the text to display to the user when they just enter '?'
///  def     : the default, ie the value to take if the user just presses return
///
{
    FILE *fp;
    SPStatus fileStat ;
    char *prmpt, *fname;
    
    size_t len = strlen(def);
    prmpt = (char *)sp_malloc(sizeof(char) * len);
    
    do {
        im_init_status(fileStat, 0) ;
        strcpy(prmpt, def);
        fname = input_string(prompt, key, help, def);
        
        if ( (fp = fopen(fname, "r")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    free(prmpt) ;
    return(fname) ;
}

char * tryWriteFile(const char *prompt, const char * key, const char * help, const char *def)
/// Function to get user input checking that the input is a valid file for writing
///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  key     : a unique bit of text (such as "AzPixelSize") which the input system will then use to store parameter values
///  help    : the text to display to the user when they just enter '?'
///  def     : the default, ie the value to take if the user just presses return
///
{
    FILE *fp;
    SPStatus fileStat ;
    char *prmpt, *fname;
    
    size_t len = strlen(def);
    prmpt = (char *)sp_malloc(sizeof(char) * len);
    
    do {
        im_init_status(fileStat, 0) ;
        strcpy(prmpt, def);
        fname = input_string(prompt, key, help, def);
        
        if ( (fp = fopen(fname, "w")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    free(prmpt) ;
    return(fname) ;
}
