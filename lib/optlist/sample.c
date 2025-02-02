/***************************************************************************
*                           OptList Usage Sample
*
*   File    : sample.cpp
*   Purpose : Demonstrates usage of optlist library.
*   Author  : Michael Dipperstein
*   Date    : July 23, 2004
*
****************************************************************************
*   HISTORY
*
*   $Id: sample.c,v 1.1.1.2 2007/09/04 04:45:42 michael Exp $
*   $Log: sample.c,v $
*   Revision 1.1.1.2  2007/09/04 04:45:42  michael
*   Added FreeOptList.
*
*   Revision 1.1.1.1  2007/08/07 05:01:48  michael
*   Initial Release
*
*
****************************************************************************
*
* Sample: A optlist library sample usage program
* Copyright (C) 2007 by Michael Dipperstein (mdipper@cs.ucsb.edu)
*
* This file is part of the optlist library.
*
* The optlist library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 3 of the
* License, or (at your option) any later version.
*
* The optlist library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
* General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
***************************************************************************/

/***************************************************************************
*                             INCLUDED FILES
***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "optlist.h"

/***************************************************************************
*                               PROTOTYPES
***************************************************************************/
char *RemovePath(char *fullPath);

/***************************************************************************
*                                FUNCTIONS
***************************************************************************/

/****************************************************************************
*   Function   : main
*   Description: This is the main function for this program, it calls
*                optlist to parse the command line input displays the
*                results of the parsing.
*   Parameters : argc - number of parameters
*                argv - parameter list
*   Effects    : parses command line parameters
*   Returned   : EXIT_SUCCESS for success, otherwise EXIT_FAILURE.
****************************************************************************/
int main(int argc, char *argv[])
{
    option_t *optList, *thisOpt;

    /* get list of command line options and their arguments */
    optList = NULL;
    optList = GetOptList(argc, argv, "a:bcd:ef?");

    /* display results of parsing */
    while (optList != NULL)
    {
        thisOpt = optList;
        optList = optList->next;

        if ('?' == thisOpt->option)
        {
            printf("Usage: %s <options>\n\n", RemovePath(argv[0]));
            printf("options:\n");
            printf("  -a : option excepting argument.\n");
            printf("  -b : option without arguments.\n");
            printf("  -c : option without arguments.\n");
            printf("  -d : option excepting argument.\n");
            printf("  -e : option without arguments.\n");
            printf("  -f : option without arguments.\n");
            printf("  -?  : print out command line options.\n\n");

            FreeOptList(thisOpt);   /* free the rest of the list */
            return EXIT_SUCCESS;
        }

        printf("found option %c\n", thisOpt->option);

        if (thisOpt->argument != NULL)
        {
            printf("\tfound argument %s", thisOpt->argument);
            printf(" at index %d\n", thisOpt->argIndex);
        }
        else
        {
            printf("\tno argument for this option\n");
        }

        free(thisOpt);    /* done with this item, free it */
    }

    return EXIT_SUCCESS;
}

/****************************************************************************
*   Function   : RemovePath
*   Description: This is function accepts a pointer to the name of a file
*                along with path information and returns a pointer to the
*                character that is not part of the path.
*   Parameters : fullPath - pointer to an array of characters containing
*                           a file name and possible path modifiers.
*   Effects    : None
*   Returned   : Returns a pointer to the first character after any path
*                information.
****************************************************************************/
char *RemovePath(char *fullPath)
{
    int i;
    char *start, *tmp;                          /* start of file name */
    const char delim[3] = {'\\', '/', ':'};     /* path deliminators */

    start = fullPath;

    /* find the first character after all file path delimiters */
    for (i = 0; i < 3; i++)
    {
        tmp = strrchr(start, delim[i]);

        if (tmp != NULL)
        {
            start = tmp + 1;
        }
    }

    return start;
}

