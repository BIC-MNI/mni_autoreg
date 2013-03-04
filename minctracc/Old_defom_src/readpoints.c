#include <stdio.h>
#include "boolean.h"

readPoints(fileName, p1, p2, numPoints)
/* this function will fill the matrix of points and put the last column
 * with 1.
 */
char *fileName; /* file name which stores the points */
float ***p1, ***p2;
int *numPoints;
/* This function will take care of reading pares of equivalent points. */
{
   FILE *fp, *fopen();
   int i;
   float **points1, **points2, **matrix();

   printf("Now reading file %s ...\n", fileName);
   if ((fp = fopen(fileName, "r")) == NULL) return(FALSE);
   /* reading a number of points in a file */
   fscanf(fp, "%d", numPoints);

   /* make approiate size of array */
   *p1 = points1 = matrix(1, *numPoints, 1, 4);
   *p2 = points2 = matrix(1, *numPoints, 1, 4);
 
   for (i=1; i<= *numPoints; i++){
      if (fscanf(fp, "%f %f %f %f %f %f", 
         &points1[i][1], &points1[i][2], &points1[i][3], 
         &points2[i][1], &points2[i][2], &points2[i][3]) == EOF) {
         printf("reading error: file name is %s", fileName);
      }
      points1[i][4] = 1; points2[i][4] = 1;
   }
   return(TRUE);
}
