#ifndef  DEF_TPS
#define  DEF_TPS

typedef struct {
  int num_points;		/* number of points used in the definition of the warp */
  int dim;			/* dimension of the warp, either ==2, or ==3           */
  double **coords;		/* matrix of coordinates in the source volume          */
				/*    that defined the warp                            */
  double **warps;		/* deformation vectors                                 */
} Thin_plate_spline;

#endif
