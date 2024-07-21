/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tool_box.h"

/****************************************************************************
*  Function     : scalloc
*  Description  : Safe calloc
*  Parameters   : number of array elements to allocate, element size, name of
                  array (for error reporting)
*  Effects      : Returns void pointer to allocated array
****************************************************************************/

void *scalloc (int n, int size, char *name)
{
  void *ptr;

  if (n <= 0) {
    fprintf (stderr, "WARNING: trying to allocate %d elements for array %s. ", n, name);
    fprintf (stderr, "returning NULL.\n");
    return NULL;
  }

  if (size <= 0) {
    fprintf (stderr, "WARNING: elements size for array %s is %d. ", name, size);
    fprintf (stderr, "returning NULL.\n");
    return NULL;
  }

  ptr = calloc (n, size);
  if (ptr == NULL) {
    fprintf (stderr, "ERROR: failed to allocate %d bytes for array %s", n * size, name);
  }

  return ptr;
}

/****************************************************************************
*  Function     : safe_free
*  Description  : Safe free memory
*  Parameters   : Pointer to allocated memory and name of structure
*  Effects      : None 
****************************************************************************/

void safe_free (void *ptr, char *name)
{
  if (ptr == NULL) {
    fprintf (stderr, "WARNING: trying to free the already NULL pointer %s!\n", name);
    return;
  }

  free (ptr);
  ptr = NULL;
}

/****************************************************************************
*  Function     : smalloc
*  Description  : Safe malloc
*  Parameters   : Number of bytes to allocate and name of allocated structure 
*  Effects      : Returns void pointer to allocated array
****************************************************************************/

void *smalloc (long n, char *name)
{
  void *ptr;

  if (n <= 0) {
    fprintf (stderr, "WARNING: trying to allocate %ld bytes for array %s. ", n, name);
    fprintf (stderr, "returning NULL.\n");
    return NULL;
  }

  ptr = malloc (n);
  if (ptr == NULL) {
    fprintf (stderr, "ERROR: failed to allocate %ld bytes for array %s", n, name);
  }

  return ptr;
}

/****************************************************************************
*  Function     : tokenize_string
*  Description  : Populates and array of string tokens from a passed string
*  Parameters   : Input string and output array of tokens
*  Effects      : Returns the number of tokens in input string and populated
                  array of tokens as second argument
****************************************************************************/

int tokenize_string (char *s, char ***tok)
{
  char test[MAX_LINE];
//  char *sep = "\t \n!=";
  char *sep = "\t \n/!=#";
  char *word;
  int count = 0;

  strncpy (test, s, MAX_LINE);

  for (word = strtok (test, sep); word; word = strtok (NULL, sep)) {
    strncpy ((*tok)[count], word, MAX_LINE);
    count++;
  }

  return count;
}
