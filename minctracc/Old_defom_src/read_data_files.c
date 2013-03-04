#include <def_standard.h>
#include <def_alloc.h>
#include <string.h>
#include <ms_iffheader.h>
#include <def_progress.h>

#include "warp_vol.h"

#include <def_voxels.h>

extern VIO_Status
   open_file(), close_file();

DATA         
   *read_iff_data();

extern Voxel_value_type 
   get_voxel_value_DATA(),
   get_voxel_value();
extern float            
   get_float_DATA();
extern void             
   put_voxel_value();

VIO_Status read_grad_data(dx,dy,dz,dxyz, name)
DATA **dx, **dy, **dz, **dxyz;
char *name;
{
   char fullname[500];
   FILE ifd;
   VIO_Status status;
   DATA *tx, *ty, *tz, *txyz;

   status = OK;

   sprintf(fullname,"%s_dx.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      tx = read_iff_data(ifd);
      close_file(ifd);
   }

   if (tx==NULL) {
      printf("problems reading in dx volume, probably not enough memory!\n ");
   } 
   else {
      *dx = tx;
   }

   sprintf(fullname,"%s_dy.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      ty = read_iff_data(ifd);
      close_file(ifd);
   }
   if (ty==NULL) {
      printf("problems reading in dx volume, probably not enough memory!\n ");
   }
   else {
      *dy = ty;
   }

   sprintf(fullname,"%s_dz.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      tz = read_iff_data(ifd);
      close_file(ifd);
   }
   if (tz==NULL) {
      printf("problems reading in dx volume, probably not enough memory!\n ");
   }
   else {
      *dz = tz;
   }

   sprintf(fullname,"%s_dxyz.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      txyz = read_iff_data(ifd);
      close_file(ifd);
   }

   if (txyz==NULL) {
      printf("problems reading in dx volume, probably not enough memory!\n ");
   }
   else {
      *dxyz = txyz;
   }

   return(status);

}

VIO_Status read_deform_data(dx,dy,dz, name)
DATA **dx, **dy, **dz;
char *name;
{
   char fullname[500];
   FILE ifd;
   VIO_Status status;
   DATA *tx, *ty, *tz;

   status = OK;

   sprintf(fullname,"%s_dx.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      tx = read_iff_data(ifd);
      close_file(ifd);
   }

   if (tx==NULL) {
      printf("problems reading in deform x volume, probably not enough memory!\n ");
   } 
   else {
      *dx = tx;
   }

   sprintf(fullname,"%s_dy.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      ty = read_iff_data(ifd);
      close_file(ifd);
   }
   if (ty==NULL) {
      printf("problems reading in deform y volume, probably not enough memory!\n ");
   }
   else {
      *dy = ty;
   }

   sprintf(fullname,"%s_dz.iff",name);
   status = open_file( fullname ,READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) {
      print_error ("input filename `%s' not found.", __FILE__, __LINE__, fullname, 0,0,0,0);
   }
   else {
      printf("%s \n",fullname);
      tz = read_iff_data(ifd);
      close_file(ifd);
   }
   if (tz==NULL) {
      printf("problems reading in deform z volume, probably not enough memory!\n ");
   }
   else {
      *dz = tz;
   }

   return(status);

}


fix_data_zeros(d)
DATA *d;
{
   int 
      row,col,slice;
   
   Voxel_value_type 
      val, val_zero;

   val_zero = get_voxel_value_DATA((float)0.0,d);

   for (slice=0; slice<d->slices; ++slice)
      for (row=0; row<d->rows; ++row)
         for (col=0; col<d->cols; ++ col) {
            val = get_voxel_value(col, row, slice, d);
            
            if (val<=1)
               put_voxel_value(val_zero, col, row, slice, d);
            
         }
}


fix_deform_zeros(dx,dy,dz)
DATA *dx, *dy, *dz;
{
   fix_data_zeros(dx);
   fix_data_zeros(dy);
   fix_data_zeros(dz);
}
