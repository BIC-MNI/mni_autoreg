#include <stdio.h>

writeTrans(filename, transMatrix)
/* This function will write the transformation matrix to the file */
char *filename;
float **transMatrix;
{
   int i, j;
   FILE *fopen(), *fp;

   fp = fopen(filename, "w");
 
   for (i=1; i<=4; i++){
      for (j=1; j<=4; j++){
         fprintf(fp, "%f ", transMatrix[i][j]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
}
