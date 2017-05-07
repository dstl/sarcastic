//
//  tryReadFile.hpp
//  sarcastic
//
//  Created by Darren Muff on 15/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef tryReadFile_hpp
#define tryReadFile_hpp

#include <stdio.h>
#include <SIlib2/SIlib2.h>
#include "colourCodes.h"

///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  help    : the text to display to the user when they just enter '?'
///  def     : the default, ie the value to take if the user just presses return
///
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def) ;

#endif /* tryReadFile_hpp */
