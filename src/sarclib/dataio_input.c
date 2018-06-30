/***************************************************************************
 *
 *       Module:    dataio_input.c
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 19/10/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions needed to get input from a user with a 
 *      prompt string and defaults.
 *
 *      Uses the GNU readline library (if available) to allow cursor editing of
 *      commands, tab completion, etc.
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  19/07/2013
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

#include "dataio.h"
#include "error_defn.h"
#include "sarclib.h"

#ifdef HAVE_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif


#define write_key(key, format, value, fname) if (fname != NULL) {                               \
    FILE * fp;                                                                                  \
                                                                                                \
    fp = fopen(fname, "a");                                                                     \
                                                                                                \
    if (fp != NULL) {                                                                           \
        if (fprintf(fp, "%s:"format"\n", key, value) < strlen(key)) {                           \
            fprintf(stderr, "Error writing key/value %s to parameter file %s\n", key, fname);   \
        }                                                                                       \
        fclose(fp);                                                                             \
    }                                                                                           \
}


static char * my_get_string(const char * prompt, const char * key, const char * help);
static char * find_key(const char * key, const char * fname);

#define N 256

// This should pull in the program name defined in image.c, and set in im_init_lib()
//
extern char * program_name;

// This will contain the input parameter file name to be used by the input_* functions
//
extern char * input_parameter_fname;

// This contains the output parameter file name to be used by the input_* functions
//
extern char * output_parameter_fname;

// Should we ask the user questions - default is we ask questions
//
extern int no_ask_questions;

static char * find_key(const char * key, const char * fname)
{
    FILE * fp;
    char * line = NULL;
    int i;
    int key_not_found = TRUE;
    
    line = calloc(N, sizeof(char));
    fp = fopen(fname, "r");
    if (fp != NULL && line != NULL)
    {
        do {
            if (fgets(line, N, fp) != NULL)
            {
                if (strncmp(line, key, strlen(key)) == 0)
                {
                    if (line[strlen(key)] == ':')
                    {
                        for(i = (int)strlen(key)+1; i < strlen(line)-1; i++)
                        {
                            line[i-(strlen(key)+1)] = line[i];
                        }
                        line[strlen(line)-1-(strlen(key)+1)] = 0;
                        key_not_found = FALSE;
                    }
                }
            }
            else
            {
                key_not_found = FALSE;
                free(line);
                line = NULL;
            }
        } while (key_not_found);
    }
    
    if (fp) fclose(fp);
    
    return (line);
}

static char * my_get_string(const char * prompt, const char * key, const char * help)
{
    char * s = NULL;
    int i;
    
    for(i = 0; i < strlen(key); i++)
    {
        if (key[i] == ':')
        {
            fprintf(stderr, "Invalid character (\":\") in key given to input_* routine\n");
            exit (7209);
        }
    }
    
    if (no_ask_questions)
    {
        printf("%s\n", prompt);
        s = calloc(5, sizeof(char));
        if (s == NULL)
        {
            fprintf(stderr, "Failed to calloc s for no ask questions!\n");
            exit(83405);
        }
        s[0] = '\n';
        s[1] = 0;
    }
    else
    {
        do {
            if (s) {
                free(s);
            }
            
#ifdef HAVE_READLINE
            s = readline(prompt);
            
            /* Read in the string, if its just a null string, then make it a newline char as
             thats what the old code expected to find if the user had just typed return.
             
             We also remove any trailing spaces (cos it gets annoying if you do a tab
             completion and then have to delete them... */
            
            if (s != NULL) {
                if (s[0] == 0) {
                    free(s);
                    s = calloc(3, sizeof(char));
                    if (s == NULL) {
                        printf("Calloc failed in my_get_string!!\n");
                        exit(78784);
                    }
                    s[0] = '\n';
                    s[1] = 0;
                } else {
                    for(i = strlen(s)-1; i > 0 && s[i] == ' '; i--) {
                        if (s[i] == ' ') {
                            s[i] = 0;
                        }
                    }
                }
            } else {
                s = calloc(3, sizeof(char));
                if (s == NULL) {
                    printf("Calloc failed in my_get_string!!\n");
                    exit(78785);
                }
                s[0] = '\n';
                s[1] = 0;
            }
            
#else
            printf ("%s", prompt);
            fflush(stdout);
            s = calloc(N, sizeof(char));
            if (s == NULL) {
                printf("Calloc failed in my_get_string!!\n");
                exit(78789);
            }
            fgets (s, N, stdin);
#endif
            
            if (strlen(s) == 1) {
                if (s[0] == '?') {
                    printf("\n%s\n\n", help);
                }
            }
        } while (strlen(s) == 1 && s[0] == '?');
    }
    
    return s;
}

// allows the user to input an int at a prompt, or select a default
//
int input_int (const char *prompt, const char * key, const char * help, int def)
{
    int val;
    char *s;
    char tprompt[N];
    char * tmp = NULL;
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        snprintf (tprompt, N, "%s [%d]: ", prompt, def);
    }
    else
    {
        snprintf (tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s != '\n')
    {
        sscanf (s, "%d", &val);
    }
    else
    {
        if (tmp == NULL)
        {
            val = def;
        }
        else
        {
            sscanf(tmp, "%d", &val);
            free(tmp);
        }
    }
    
    free(s);
    
    write_key(key, "%d", val, output_parameter_fname);
    
    return (val);
}

// allows the user to input an int at a prompt, or select a default
//
int64_t input_int64 (const char *prompt, const char * key, const char * help, int64_t def)
{
    int64_t val;
    char *s;
    char tprompt[N];
    char * tmp = NULL;
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        snprintf (tprompt, N, "%s [%"PRId64"]: ", prompt, def);
    }
    else
    {
        snprintf (tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s != '\n')
    {
        sscanf (s, "%"PRId64, &val);
    }
    else
    {
        if (tmp == NULL)
        {
            val = def;
        }
        else
        {
            sscanf(tmp, "%"PRId64, &val);
            free(tmp);
        }
    }
    
    free(s);
    
    write_key(key, "%"PRId64, val, output_parameter_fname);
    
    return (val);
}

// allows the user to input a yes (=1)/no(=0) at a prompt, or select a default
//
int input_yesno (const char *prompt, const char * key, const char * help, int def)
{
    int val;
    char * s;
    char def_prompt;
    char tprompt[N];
    char * tmp = NULL;
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        if (def) {
            def_prompt = 'Y';
        } else {
            def_prompt = 'N';
        }
        snprintf(tprompt, N, "%s [%c]: ", prompt, def_prompt);
    }
    else
    {
        snprintf(tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s != '\n')
    {
        val = (*s == 'y' || *s == 'Y');
    }
    else
    {
        if (tmp == NULL)
        {
            val = def;
        }
        else
        {
            val = (*tmp == 'y' || *tmp == 'Y');
            free(tmp);
        }
    }
    
    free(s);
    
    write_key(key, "%c", ((val) ? 'Y' : 'N'), output_parameter_fname);
    
    return (val);
}

// allows the user to input an double at a prompt, or select a default
//
double input_dbl (const char *prompt, const char * key, const char * help, double def)
{
    double val;
    char *s;
    char * tmp = NULL;
    char tprompt[N];
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        snprintf (tprompt, N, "%s [%.10g]: ", prompt, def);
    }
    else
    {
        snprintf (tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s != '\n')
    {
        sscanf (s, "%lf", &val);
    }
    else
    {
        if (tmp == NULL)
        {
            val = def;
        }
        else
        {
            sscanf(tmp, "%lf", &val);
            free(tmp);
        }
    }
    
    free(s);
    
    write_key(key, "%.14g", val, output_parameter_fname);
    
    return (val);
}

// allows the user to input an SPVector at a prompt, or select a default
//
SPVector input_vect (const char *prompt, const char * key, const char * help, SPVector def)      
{
    SPVector val;
    char * s;
    char * tmp = NULL;
    char tprompt[N];
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        snprintf (tprompt, N, "%s [%.8g %.8g %.8g]: ", prompt, def.x, def.y, def.z);
    }
    else
    {
        snprintf (tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s != '\n')
    {
        sscanf (s, "%lf %lf %lf", &val.x,  &val.y,  &val.z);
    }
    else
    {
        if (tmp == NULL)
        {
            val = def;
        }
        else
        {
            sscanf (tmp, "%lf %lf %lf", &val.x,  &val.y,  &val.z);
            free(tmp);
        }
    }
    
    tmp = calloc(3 *30, sizeof(char));
    snprintf(tmp, 3*30, "%.14g %.14g %.14g", val.x, val.y, val.z);
    write_key(key, "%s", tmp, output_parameter_fname);
    free(tmp);
    
    return (val);
}

// allows the user to input a string at a prompt, or select a default
//
char * input_string (const char *prompt, const char * key, const char * help, const char *def)
{
    char * s;
    char tprompt[N];
    char * tmp = NULL;
    
    if (input_parameter_fname != NULL)
    {
        tmp = find_key(key, input_parameter_fname);
    }
    
    if (tmp == NULL)
    {
        snprintf(tprompt, N, "%s [%s]: ", prompt, def);
    }
    else
    {
        snprintf (tprompt, N, "%s [%s]: ", prompt, tmp);
    }
    
    s = my_get_string(tprompt, key, help);
    
    if (*s == '\n') 
    {
        free(s);
        if (tmp == NULL)
        {
            s = strdup(def);
        }
        else
        {
            s = tmp;
        }
    }
    
    
    if (s[strlen(s)-1] == '\n') {
        s[strlen(s)-1] = 0;
    }
    
    write_key(key, "%s", s, output_parameter_fname);
    
    return (s);
}

