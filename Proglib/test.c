#include <stdio.h>
#include <string.h>
#include <minc.h>
#include <minc_def.h>
#include <volume_io.h>

#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif


/* Main program */

int get_file_dim_names(char *filename, int *dim, char ***dim_names)
{
   int mincid, varid, dimid;
   int ndims, dims[MAX_VAR_DIMS];
   int idim;

   char name[MAX_NC_NAME];


   ncopts = NC_VERBOSE | NC_FATAL;

   /* Open the file */
   mincid = ncopen(filename, NC_NOWRITE);


   if (ncinquire(mincid, &ndims, NULL, NULL, NULL) == MI_ERROR) {
     printf ("Error\n"); exit(FALSE);
   }
   else {
     ALLOC(*dim_names,ndims);
     *dim = ndims;
     for (idim=0; idim<ndims; idim++) {
       if (ncdiminq(mincid, idim, name, NULL) == MI_ERROR) {
	 printf ("Error\n"); exit(FALSE);
       }
       else {
	 (void) printf("%s ", name);
	 ALLOC( (*dim_names)[idim],(strlen(name)+3));
	 strcat((*dim_names)[idim],"MI");
	 strcat((*dim_names)[idim],name);
	 
       }
     }
     (void) printf("\n");
   }
   
   /* Close the file */
   (void) ncclose(mincid);
   

}

int main(int argc, char *argv[]) {
 char **names;
 int i,dim;

 get_file_dim_names(argv[1],&dim,&names);

 printf("%d\n",dim);

 for(i=0; i<dim; i++) {
   printf ("%d: %s\n",i,names[i]);
 }

 FREE2D(names);


}
