#define DO_TRANSFORM(result, transformation, coord) \
   general_transform_point(transformation, \
      Point_x(coord), Point_y(coord), Point_z(coord), \
      &(Point_x(result)), &(Point_y(result)), &(Point_z(result)) )

#define IS_LINEAR(transformation) \
   (get_transform_type(transformation)==LINEAR)


#define VOLUME_VALUE(volume, ind0, ind1, ind2, value) \
{ \
   long _offset_; \
 \
   _offset_ = ((ind0)*volume->size[1] + (ind1))*volume->size[2] + (ind2); \
   switch (volume->datatype) { \
   case NC_BYTE: \
      if (volume->is_signed) \
         value = *((signed char *) volume->data + _offset_); \
      else \
         value = *((unsigned char *) volume->data + _offset_); \
      break; \
   case NC_SHORT: \
      if (volume->is_signed) \
         value = *((signed short *) volume->data + _offset_); \
      else \
         value = *((unsigned short *) volume->data + _offset_); \
      break; \
   case NC_LONG: \
      if (volume->is_signed) \
         value = *((signed long *) volume->data + _offset_); \
      else \
         value = *((unsigned long *) volume->data + _offset_); \
      break; \
   case NC_FLOAT: \
      value = *((float *) volume->data + _offset_); \
      break; \
   case NC_DOUBLE: \
      value = *((double *) volume->data + _offset_); \
      break; \
   } \
}

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   (*(main_args.interpolant)) (volume, coord, result)

#ifndef DEBUG_PRINT
#   define DEBUG_PRINT(str) if (main_args.flags.debug) (void) fprintf (stderr,  str  );
#   define DEBUG_PRINT1(str,a1) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1 );
#   define DEBUG_PRINT2(str,a1,a2) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2 );
#   define DEBUG_PRINT3(str,a1,a2,a3) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3 );
#   define DEBUG_PRINT4(str,a1,a2,a3,a4) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4 );
#   define DEBUG_PRINT5(str,a1,a2,a3,a4,a5) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5 );
#   define DEBUG_PRINT6(str,a1,a2,a3,a4,a5,a6) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6 );
#   define DEBUG_PRINT7(str,a1,a2,a3,a4,a5,a6,a7) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6,a7 );
#endif

#ifdef INTERPOLATE
#  undef INTERPOLATE
#endif

#define INTERPOLATE(volume, coord, result) \
   (*volume->interpolant) (volume, coord, result)

/* These are needed if linking against (eg.) the GNU C Library, since
 * it has only the ANSI-named functions (sqrtf, etc.) and not the 
 * more traditional fsqrt, etc.  Added by GPW 95/2/9.
 */

#ifdef STRICT_ANSI_LIB
# define fsqrt sqrtf
# define fcos cosf
# define fsin sinf
# define facos acosf
# define fasin asinf
#endif


