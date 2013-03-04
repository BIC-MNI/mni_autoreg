/* ----------------------------- MNI Header -----------------------------------
@NAME       : quaternion.h
@DESCRIPTION: Header file for quaternion.c
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
---------------------------------------------------------------------------- */


/*  ------------------------ Function prototypes  ------------------------ */
void vcopy(double *copy, double *v);


void vadd(double *add, double *v1, double *v2);


void vscale(double *v, double scale);


void vcross(double *cross, double *v1, double *v2);


double vdot(double *v1, double *v2);


double vlength(double *v);


double veuc(double *p1, double *p2);

void vnormal(double *v);


void axis_to_quat(double vec[3], double phi, double quat[4]);


void quat_to_axis(double vec[3], double *phi, double quat[4]);


void add_quats(double q1[4], double q2[4], double dest[4]);

void build_rotmatrix(float **m, double *quat);

void extract_quaternions(float **m, double *quat);



