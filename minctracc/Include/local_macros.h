
#ifndef LOCAL_MACROS_H
#define LOCAL_MACROS_H

/* a bunch of (evil) local macros */


/* from volume_io.h */
#define  ALLOC( ptr, n_items )                                                \
               ptr = alloc_memory_1d( (size_t) (n_items),                     \
                         sizeof(*(ptr)) _ALLOC_SOURCE_LINE )

#define  REALLOC( ptr, n_items )                                              \
           realloc_memory( (void **) &(ptr), (size_t) (n_items),              \
                         sizeof(*(ptr)) _ALLOC_SOURCE_LINE )


#define  FREE( ptr )                                                          \
   free_memory_1d( (void **) &(ptr) _ALLOC_SOURCE_LINE )

#define  ALLOC2D( ptr, n1, n2 )                                               \
               ptr = alloc_memory_2d( (size_t) (n1), (size_t) (n2),           \
                          sizeof(**(ptr)) _ALLOC_SOURCE_LINE )

#define  FREE2D( ptr )                                                        \
         free_memory_2d( (void ***) &(ptr) _ALLOC_SOURCE_LINE )


#define DO_TRANSFORM(result, transformation, coord) \
   general_transform_point(transformation, \
      Point_x(coord), Point_y(coord), Point_z(coord), \
      &(Point_x(result)), &(Point_y(result)), &(Point_z(result)) )

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   (*(main_args.interpolant)) (volume, coord, result)

#ifndef DEBUG_PRINT
#   define DEBUG_PRINT(str) if (main_args.flags.debug) (void) fprintf (stderr,  str  );
#   define DEBUG_PRINT1(str,a1) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1 );
#   define DEBUG_PRINT2(str,a1,a2) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2 );
#   define DEBUG_PRINT3(str,a1,a2,a3) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3 );
#   define DEBUG_PRINT4(str,a1,a2,a3,a4) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4 );
#endif


#endif /* LOCAL_MACROS_H */
