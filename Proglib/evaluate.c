#include  <volume_io.h>
#include  "splines.h"


public  void   evaluate_3D_volume_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz );

public  void   special_evaluate_3D_volume_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    BOOLEAN        x_use_higher,
    BOOLEAN        y_use_higher,
    BOOLEAN        z_use_higher,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z );

public  void   evaluate_3D_slice_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_yy );

public  void   evaluate_3D_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz );

public  void   evaluate_3D_slice(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_yy );


private  void   trilinear_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    BOOLEAN        x_use_higher,
    BOOLEAN        y_use_higher,
    BOOLEAN        z_use_higher,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_3D_volume_in_world
@INPUT      : volume
              x
              y
              z
              degrees_continuity - 0 = linear, 2 = cubic
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_z
@RETURNS    : 
@DESCRIPTION: Takes a world space position and evaluates the value within
              the volume.
              If deriv_x is not a null pointer, then the 3 derivatives are
              passed back.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

public  void   evaluate_3D_volume_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz )
{
    Real      ignore, dxx, dxy, dxz, dyy, dyz, dzz;
    Real      voxel[MAX_DIMENSIONS];
    Real      txx, txy, txz;
    Real      tyx, tyy, tyz;
    Real      tzx, tzy, tzz;

    if( get_volume_n_dimensions(volume) != 3 )
    {
        handle_internal_error(
                 "evaluate_3D_volume_in_world: volume must be 3D.\n" );
    }

    convert_world_to_voxel( volume, x, y, z, voxel );

    evaluate_3D_volume( volume, voxel[X], voxel[Y], voxel[Z],
                        degrees_continuity, value, deriv_x, deriv_y, deriv_z,
                        deriv_xx, deriv_xy, deriv_xz,
                        deriv_yy, deriv_yz, deriv_zz );

    if( deriv_x != (Real *) 0 )
    {
        convert_voxel_normal_vector_to_world( volume,
                                              *deriv_x, *deriv_y, *deriv_z,
                                              deriv_x, deriv_y, deriv_z );
    }

    if( deriv_xx != (Real *) 0 )
    {
        dxx = *deriv_xx;
        dxy = *deriv_xy;
        dxz = *deriv_xz;
        dyy = *deriv_yy;
        dyz = *deriv_yz;
        dzz = *deriv_zz;
        convert_voxel_normal_vector_to_world( volume,
                                              dxx, dxy, dxz,
                                              &txx, &txy, &txz );
        convert_voxel_normal_vector_to_world( volume,
                                              dxy, dyy, dyz,
                                              &tyx, &tyy, &tyz );
        convert_voxel_normal_vector_to_world( volume,
                                              dxz, dyz, dzz,
                                              &tzx, &tzy, &tzz );

        convert_voxel_normal_vector_to_world( volume,
                                              txx, tyx, tzx,
                                              deriv_xx, &ignore, &ignore );
        convert_voxel_normal_vector_to_world( volume,
                                              txy, tyy, tzy,
                                              deriv_xy, deriv_yy, &ignore );
        convert_voxel_normal_vector_to_world( volume,
                                              txz, tyz, tzz,
                                              deriv_xz, deriv_yz, deriv_zz );
    }
}

public  void   special_evaluate_3D_volume_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    BOOLEAN        x_use_higher,
    BOOLEAN        y_use_higher,
    BOOLEAN        z_use_higher,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z )
{
    Real      voxel[MAX_DIMENSIONS];
    int       sizes[MAX_DIMENSIONS];

    convert_world_to_voxel( volume, x, y, z, voxel );

    get_volume_sizes( volume, sizes );
    if( x < -0.5 || x > (Real) sizes[X] - 0.5 ||
        y < -0.5 || y > (Real) sizes[Y] - 0.5 ||
        z < -0.5 || z > (Real) sizes[Z] - 0.5 )
    {
        *value = get_volume_voxel_min( volume );
        if( deriv_x != (Real *) NULL )
        {
            *deriv_x = 0.0;
            *deriv_y = 0.0;
            *deriv_z = 0.0;
        }
        return;
    }

    trilinear_interpolate_volume( volume, voxel[X], voxel[Y], voxel[Z],
                                  x_use_higher, y_use_higher, z_use_higher,
                                  value, deriv_x, deriv_y, deriv_z );

    if( deriv_x != (Real *) 0 )
    {
        convert_voxel_normal_vector_to_world( volume,
                                              *deriv_x, *deriv_y, *deriv_z,
                                              deriv_x, deriv_y, deriv_z );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_3D_slice_in_world
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_xx
              deriv_xy
              deriv_yy
@RETURNS    : 
@DESCRIPTION: Takes a world space position and evaluates the value within
              the volume by bilinear interpolation within the slice.
              If deriv_x is not a null pointer, then the derivatives are passed
              back.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

public  void   evaluate_3D_slice_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_yy )
{
    Real      ignore, dxx, dxy, dyy;
    Real      voxel[MAX_DIMENSIONS];
    Real      txx, txy, txz;
    Real      tyx, tyy, tyz;
    Real      tzx, tzy, tzz;

    if( get_volume_n_dimensions(volume) != 3 )
    {
        handle_internal_error(
                 "evaluate_3D_slice_in_world: volume must be 3D.\n" );
    }

    convert_world_to_voxel( volume, x, y, z, voxel );

    evaluate_3D_slice( volume, voxel[X], voxel[Y], voxel[Z],
                       value, deriv_x, deriv_y,
                       deriv_xx, deriv_xy, deriv_yy );

    if( deriv_x != (Real *) 0 )
    {
        convert_voxel_normal_vector_to_world( volume,
                                              *deriv_x, *deriv_y, 0.0,
                                              deriv_x, deriv_y, &ignore );
    }

    if( deriv_xx != (Real *) 0 )
    {
        dxx = *deriv_xx;
        dxy = *deriv_xy;
        dyy = *deriv_yy;
        convert_voxel_normal_vector_to_world( volume,
                                              dxx, dxy, 0.0,
                                              &txx, &txy, &txz );
        convert_voxel_normal_vector_to_world( volume,
                                              dxy, dyy, 0.0,
                                              &tyx, &tyy, &tyz );
        convert_voxel_normal_vector_to_world( volume,
                                              0.0, 0.0, 0.0,
                                              &tzx, &tzy, &tzz );

        convert_voxel_normal_vector_to_world( volume,
                                              txx, tyx, tzx,
                                              deriv_xx, &ignore, &ignore );
        convert_voxel_normal_vector_to_world( volume,
                                              txy, tyy, tzy,
                                              deriv_xy, deriv_yy, &ignore );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : triconstant_interpolate_volume
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_z
@RETURNS    : 
@DESCRIPTION: Returns the value within the volume, assuming constant voxels,
              (step function).  Derivative is approximated by neighbours.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  void   triconstant_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z )
{
    int                i, j, k, sizes[MAX_DIMENSIONS];
    Real               prev, next, dx, dy, dz, voxel_value;

    get_volume_sizes( volume, sizes );

    i = ROUND( x );
    if( i == sizes[X] )
        i = sizes[X]-1;
    j = ROUND( y );
    if( j == sizes[Y] )
        j = sizes[Y]-1;
    k = ROUND( z );
    if( k == sizes[Z] )
        k = sizes[Z]-1;

    GET_VALUE_3D( voxel_value, volume, i, j, k );

    if( value != (Real *) 0 )
        *value = voxel_value;

    if( deriv_x != (Real *) NULL )
    {
        /* --- get derivative wrt x */

        dx = 0;
        if( i == 0 )
            prev = voxel_value;
        else
        {
            GET_VALUE_3D( prev, volume, i-1, j, k );
            ++dx;
        }

        if( i == sizes[X]-1 )
            next = voxel_value;
        else
        {
            GET_VALUE_3D( next, volume, i+1, j, k );
            ++dx;
        }

        if( dx == 0 )
            *deriv_x = 0.0;
        else
            *deriv_x = (next - prev) / (Real) dx;

        /* --- get derivative wrt y */

        dy = 0;
        if( j == 0 )
            prev = voxel_value;
        else
        {
            GET_VALUE_3D( prev, volume, i, j-1, k );
            ++dy;
        }

        if( j == sizes[Y]-1 )
            next = voxel_value;
        else
        {
            GET_VALUE_3D( next, volume, i, j+1, k );
            ++dy;
        }

        if( dy == 0 )
            *deriv_y = 0.0;
        else
            *deriv_y = (next - prev) / (Real) dy;

        /* --- get derivative wrt z */

        dz = 0;
        if( k == 0 )
            prev = voxel_value;
        else
        {
            GET_VALUE_3D( prev, volume, i, j, k-1 );
            ++dz;
        }

        if( k == sizes[Z]-1 )
            next = voxel_value;
        else
        {
            GET_VALUE_3D( next, volume, i, j, k+1 );
            ++dz;
        }

        if( dz == 0 )
            *deriv_z = 0.0;
        else
            *deriv_z = (next - prev) / (Real) dz;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : trilinear_interpolate_volume
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_z
@RETURNS    : 
@DESCRIPTION: Returns the value within the volume, assuming trilinear
              interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  void   trilinear_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    BOOLEAN        x_use_higher,
    BOOLEAN        y_use_higher,
    BOOLEAN        z_use_higher,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z )
{
    int                i, j, k, sizes[MAX_DIMENSIONS];
    Real               u, v, w;
    Real               c000, c001, c010, c011, c100, c101, c110, c111;
    Real               c00, c01, c10, c11;
    Real               c0, c1;
    Real               du00, du01, du10, du11, du0, du1;
    Real               dv0, dv1;

    get_volume_sizes( volume, sizes );

    if( x == (Real) (sizes[X]-1) )
    {
        i = sizes[X]-2;
        u = 1.0;
    }
    else
    {
        i = (int) x;
        u = FRACTION( x );
        if( !x_use_higher && u == 0.0 && i > 0 )
        {
            --i;
            u = 1.0;
        }
    }

    if( y == (Real) (sizes[Y]-1) )
    {
        j = sizes[Y]-2;
        v = 1.0;
    }
    else
    {
        j = (int) y;
        v = FRACTION( y );
        if( !y_use_higher && v == 0.0 && j > 0 )
        {
            --j;
            v = 1.0;
        }
    }

    if( z == (Real) (sizes[Z]-1) )
    {
        k = sizes[Z]-2;
        w = 1.0;
    }
    else
    {
        k = (int) z;
        w = FRACTION( z );
        if( !z_use_higher && w == 0.0 && k > 0 )
        {
            --k;
            w = 1.0;
        }
    }

    GET_VALUE_3D( c000, volume, i,   j,   k );
    GET_VALUE_3D( c001, volume, i,   j,   k+1 );
    GET_VALUE_3D( c010, volume, i,   j+1, k );
    GET_VALUE_3D( c011, volume, i,   j+1, k+1 );
    GET_VALUE_3D( c100, volume, i+1, j,   k );
    GET_VALUE_3D( c101, volume, i+1, j,   k+1 );
    GET_VALUE_3D( c110, volume, i+1, j+1, k );
    GET_VALUE_3D( c111, volume, i+1, j+1, k+1 );

    du00 = c100 - c000;
    du01 = c101 - c001;
    du10 = c110 - c010;
    du11 = c111 - c011;

    c00 = c000 + u * du00;
    c01 = c001 + u * du01;
    c10 = c010 + u * du10;
    c11 = c011 + u * du11;

    dv0 = c10 - c00;
    dv1 = c11 - c01;

    c0 = c00 + v * dv0;
    c1 = c01 + v * dv1;

    if( value != (Real *) 0 )
        *value = INTERPOLATE( w, c0, c1 );

    if( deriv_x != (Real *) 0 )
    {
        du0 = INTERPOLATE( v, du00, du10 );
        du1 = INTERPOLATE( v, du01, du11 );

        *deriv_x = INTERPOLATE( w, du0, du1 );
        *deriv_y = INTERPOLATE( w, dv0, dv1 );
        *deriv_z = (c1 - c0);
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : triquadratic_interpolate_volume
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_z
@RETURNS    : 
@DESCRIPTION: Returns the value within the volume, assuming trilinear
              interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  void   triquadratic_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz )
{
    int      i, j, k;
    Real     tx, ty, tz, u, v, w, dummy;
    Real     c000, c001, c002, c010, c011, c012, c020, c021, c022;
    Real     c100, c101, c102, c110, c111, c112, c120, c121, c122;
    Real     c200, c201, c202, c210, c211, c212, c220, c221, c222;
    int      sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    tx = x + 0.5;
    ty = y + 0.5;
    tz = z + 0.5;

    if( x == (Real) sizes[X] - 1.5 )
    {
        i = sizes[X]-2;
        u = 1.0;
    }
    else
    {
        i = (int) tx;
        u = FRACTION( tx );
    }

    if( y == (Real) sizes[Y] - 1.5 )
    {
        j = sizes[Y]-2;
        v = 1.0;
    }
    else
    {
        j = (int) ty;
        v = FRACTION( ty );
    }

    if( z == (Real) sizes[Z] - 1.5 )
    {
        k = sizes[Z]-2;
        w = 1.0;
    }
    else
    {
        k = (int) tz;
        w = FRACTION( tz );
    }

    GET_VALUE_3D( c000, volume, i-1, j-1, k-1 );
    GET_VALUE_3D( c001, volume, i-1, j-1, k+0 );
    GET_VALUE_3D( c002, volume, i-1, j-1, k+1 );
    GET_VALUE_3D( c010, volume, i-1, j+0, k-1 );
    GET_VALUE_3D( c011, volume, i-1, j+0, k+0 );
    GET_VALUE_3D( c012, volume, i-1, j+0, k+1 );
    GET_VALUE_3D( c020, volume, i-1, j+1, k-1 );
    GET_VALUE_3D( c021, volume, i-1, j+1, k+0 );
    GET_VALUE_3D( c022, volume, i-1, j+1, k+1 );

    GET_VALUE_3D( c100, volume, i+0, j-1, k-1 );
    GET_VALUE_3D( c101, volume, i+0, j-1, k+0 );
    GET_VALUE_3D( c102, volume, i+0, j-1, k+1 );
    GET_VALUE_3D( c110, volume, i+0, j+0, k-1 );
    GET_VALUE_3D( c111, volume, i+0, j+0, k+0 );
    GET_VALUE_3D( c112, volume, i+0, j+0, k+1 );
    GET_VALUE_3D( c120, volume, i+0, j+1, k-1 );
    GET_VALUE_3D( c121, volume, i+0, j+1, k+0 );
    GET_VALUE_3D( c122, volume, i+0, j+1, k+1 );

    GET_VALUE_3D( c200, volume, i+1, j-1, k-1 );
    GET_VALUE_3D( c201, volume, i+1, j-1, k+0 );
    GET_VALUE_3D( c202, volume, i+1, j-1, k+1 );
    GET_VALUE_3D( c210, volume, i+1, j+0, k-1 );
    GET_VALUE_3D( c211, volume, i+1, j+0, k+0 );
    GET_VALUE_3D( c212, volume, i+1, j+0, k+1 );
    GET_VALUE_3D( c220, volume, i+1, j+1, k-1 );
    GET_VALUE_3D( c221, volume, i+1, j+1, k+0 );
    GET_VALUE_3D( c222, volume, i+1, j+1, k+1 );

    if( value != (Real *) 0 )
    {
        QUADRATIC_TRIVAR(c, u, v, w, *value );
    }

    if( deriv_x != (Real *) 0 )
    {
        QUADRATIC_TRIVAR_DERIV(c, u, v, w, dummy, *deriv_x, *deriv_y, *deriv_z );
    }

    if( deriv_xx != (Real *) 0 )
    {
        QUADRATIC_TRIVAR_DERIV2(c, u, v, w,
                                 dummy, dummy, dummy, dummy,
                                 *deriv_xx, *deriv_xy, *deriv_xz,
                                 *deriv_yy, *deriv_yz, *deriv_zz );
    }

#ifdef  lint
    if( dummy == 0.0 ) {}
#endif
}

private  void   tricubic_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz )
{
    int                i, j, k;
    Real               u, v, w, dummy;
    Real               c000, c001, c002, c003, c010, c011, c012, c013;
    Real               c020, c021, c022, c023, c030, c031, c032, c033;
    Real               c100, c101, c102, c103, c110, c111, c112, c113;
    Real               c120, c121, c122, c123, c130, c131, c132, c133;
    Real               c200, c201, c202, c203, c210, c211, c212, c213;
    Real               c220, c221, c222, c223, c230, c231, c232, c233;
    Real               c300, c301, c302, c303, c310, c311, c312, c313;
    Real               c320, c321, c322, c323, c330, c331, c332, c333;
    int                sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    if( x == (Real) sizes[X] - 1.5 )
    {
        i = sizes[X]-2;
        u = 1.0;
    }
    else
    {
        i = (int) x;
        u = FRACTION( x );
    }

    if( y == (Real) sizes[Y] - 1.5 )
    {
        j = sizes[Y]-2;
        v = 1.0;
    }
    else
    {
        j = (int) y;
        v = FRACTION( y );
    }

    if( z == (Real) sizes[Z] - 1.5 )
    {
        k = sizes[Z]-2;
        w = 1.0;
    }
    else
    {
        k = (int) z;
        w = FRACTION( z );
    }

    GET_VALUE_3D( c000, volume, i-1, j-1, k-1 );
    GET_VALUE_3D( c001, volume, i-1, j-1, k+0 );
    GET_VALUE_3D( c002, volume, i-1, j-1, k+1 );
    GET_VALUE_3D( c003, volume, i-1, j-1, k+2 );
    GET_VALUE_3D( c010, volume, i-1, j+0, k-1 );
    GET_VALUE_3D( c011, volume, i-1, j+0, k+0 );
    GET_VALUE_3D( c012, volume, i-1, j+0, k+1 );
    GET_VALUE_3D( c013, volume, i-1, j+0, k+2 );
    GET_VALUE_3D( c020, volume, i-1, j+1, k-1 );
    GET_VALUE_3D( c021, volume, i-1, j+1, k+0 );
    GET_VALUE_3D( c022, volume, i-1, j+1, k+1 );
    GET_VALUE_3D( c023, volume, i-1, j+1, k+2 );
    GET_VALUE_3D( c030, volume, i-1, j+2, k-1 );
    GET_VALUE_3D( c031, volume, i-1, j+2, k+0 );
    GET_VALUE_3D( c032, volume, i-1, j+2, k+1 );
    GET_VALUE_3D( c033, volume, i-1, j+2, k+2 );

    GET_VALUE_3D( c100, volume, i+0, j-1, k-1 );
    GET_VALUE_3D( c101, volume, i+0, j-1, k+0 );
    GET_VALUE_3D( c102, volume, i+0, j-1, k+1 );
    GET_VALUE_3D( c103, volume, i+0, j-1, k+2 );
    GET_VALUE_3D( c110, volume, i+0, j+0, k-1 );
    GET_VALUE_3D( c111, volume, i+0, j+0, k+0 );
    GET_VALUE_3D( c112, volume, i+0, j+0, k+1 );
    GET_VALUE_3D( c113, volume, i+0, j+0, k+2 );
    GET_VALUE_3D( c120, volume, i+0, j+1, k-1 );
    GET_VALUE_3D( c121, volume, i+0, j+1, k+0 );
    GET_VALUE_3D( c122, volume, i+0, j+1, k+1 );
    GET_VALUE_3D( c123, volume, i+0, j+1, k+2 );
    GET_VALUE_3D( c130, volume, i+0, j+2, k-1 );
    GET_VALUE_3D( c131, volume, i+0, j+2, k+0 );
    GET_VALUE_3D( c132, volume, i+0, j+2, k+1 );
    GET_VALUE_3D( c133, volume, i+0, j+2, k+2 );

    GET_VALUE_3D( c200, volume, i+1, j-1, k-1 );
    GET_VALUE_3D( c201, volume, i+1, j-1, k+0 );
    GET_VALUE_3D( c202, volume, i+1, j-1, k+1 );
    GET_VALUE_3D( c203, volume, i+1, j-1, k+2 );
    GET_VALUE_3D( c210, volume, i+1, j+0, k-1 );
    GET_VALUE_3D( c211, volume, i+1, j+0, k+0 );
    GET_VALUE_3D( c212, volume, i+1, j+0, k+1 );
    GET_VALUE_3D( c213, volume, i+1, j+0, k+2 );
    GET_VALUE_3D( c220, volume, i+1, j+1, k-1 );
    GET_VALUE_3D( c221, volume, i+1, j+1, k+0 );
    GET_VALUE_3D( c222, volume, i+1, j+1, k+1 );
    GET_VALUE_3D( c223, volume, i+1, j+1, k+2 );
    GET_VALUE_3D( c230, volume, i+1, j+2, k-1 );
    GET_VALUE_3D( c231, volume, i+1, j+2, k+0 );
    GET_VALUE_3D( c232, volume, i+1, j+2, k+1 );
    GET_VALUE_3D( c233, volume, i+1, j+2, k+2 );

    GET_VALUE_3D( c300, volume, i+2, j-1, k-1 );
    GET_VALUE_3D( c301, volume, i+2, j-1, k+0 );
    GET_VALUE_3D( c302, volume, i+2, j-1, k+1 );
    GET_VALUE_3D( c303, volume, i+2, j-1, k+2 );
    GET_VALUE_3D( c310, volume, i+2, j+0, k-1 );
    GET_VALUE_3D( c311, volume, i+2, j+0, k+0 );
    GET_VALUE_3D( c312, volume, i+2, j+0, k+1 );
    GET_VALUE_3D( c313, volume, i+2, j+0, k+2 );
    GET_VALUE_3D( c320, volume, i+2, j+1, k-1 );
    GET_VALUE_3D( c321, volume, i+2, j+1, k+0 );
    GET_VALUE_3D( c322, volume, i+2, j+1, k+1 );
    GET_VALUE_3D( c323, volume, i+2, j+1, k+2 );
    GET_VALUE_3D( c330, volume, i+2, j+2, k-1 );
    GET_VALUE_3D( c331, volume, i+2, j+2, k+0 );
    GET_VALUE_3D( c332, volume, i+2, j+2, k+1 );
    GET_VALUE_3D( c333, volume, i+2, j+2, k+2 );

    if( deriv_xx != (Real *) 0 )
    {
        CUBIC_TRIVAR_DERIV2(c, u, v, w, *value,
                             *deriv_x, *deriv_y, *deriv_z,
                             *deriv_xx, *deriv_xy, *deriv_xz,
                             *deriv_yy, *deriv_yz, *deriv_zz );
    }
    else
    {
        if( value != (Real *) 0 )
        {
            CUBIC_TRIVAR(c, u, v, w, *value );
        }

        if( deriv_x != (Real *) 0 )
        {
            CUBIC_TRIVAR_DERIV(c, u, v, w,
                                dummy, *deriv_x, *deriv_y, *deriv_z );
        }
    }

#ifdef  lint
    if( dummy == 0.0 ) {}
#endif
}

private  void   bicubic_interpolate_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_yy )
{
    int                i, j, k;
    Real               u, v, dummy;
    Real               c00, c01, c02, c03, c10, c11, c12, c13;
    Real               c20, c21, c22, c23, c30, c31, c32, c33;
    int                sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    if( x == (Real) sizes[X] - 1.5 )
    {
        i = sizes[X]-2;
        u = 1.0;
    }
    else
    {
        i = (int) x;
        u = FRACTION( x );
    }

    if( y == (Real) sizes[Y] - 1.5 )
    {
        j = sizes[Y]-2;
        v = 1.0;
    }
    else
    {
        j = (int) y;
        v = FRACTION( y );
    }

    if( z == (Real) sizes[Z] - 1.5 )
    {
        k = sizes[Z]-2;
    }
    else
    {
        k = (int) z;
    }

    GET_VALUE_3D( c00, volume, i-1, j-1, k );
    GET_VALUE_3D( c01, volume, i-1, j+0, k );
    GET_VALUE_3D( c02, volume, i-1, j+1, k );
    GET_VALUE_3D( c03, volume, i-1, j+2, k );
    GET_VALUE_3D( c10, volume, i+0, j-1, k );
    GET_VALUE_3D( c11, volume, i+0, j+0, k );
    GET_VALUE_3D( c12, volume, i+0, j+1, k );
    GET_VALUE_3D( c13, volume, i+0, j+2, k );
    GET_VALUE_3D( c20, volume, i+1, j-1, k );
    GET_VALUE_3D( c21, volume, i+1, j+0, k );
    GET_VALUE_3D( c22, volume, i+1, j+1, k );
    GET_VALUE_3D( c23, volume, i+1, j+2, k );
    GET_VALUE_3D( c30, volume, i+2, j-1, k );
    GET_VALUE_3D( c31, volume, i+2, j+0, k );
    GET_VALUE_3D( c32, volume, i+2, j+1, k );
    GET_VALUE_3D( c33, volume, i+2, j+2, k );

    if( deriv_xx != (Real *) 0 )
    {
        CUBIC_BIVAR_DERIV2(c, u, v, *value, *deriv_x, *deriv_y,
                            *deriv_xx, *deriv_xy, *deriv_yy );
    }
    else
    {
        if( value != (Real *) 0 )
        {
            CUBIC_BIVAR(c, u, v, *value );
        }

        if( deriv_x != (Real *) 0 )
        {
            CUBIC_BIVAR_DERIV(c, u, v, dummy, *deriv_x, *deriv_y );
        }
    }

#ifdef  lint
    if( dummy == 0.0 ) {}
#endif
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_3D_volume
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_z
@RETURNS    : 
@DESCRIPTION: Takes a voxel space position and evaluates the value within
              the volume by trilinear interpolation.
              If deriv_x is not a null pointer, then the 3 derivatives are passed
              back.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

public  void   evaluate_3D_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_z,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_xz,
    Real           *deriv_yy,
    Real           *deriv_yz,
    Real           *deriv_zz )
{
    int         sizes[MAX_DIMENSIONS];

    if( get_volume_n_dimensions(volume) != 3 )
    {
        handle_internal_error( "evaluate_3D_volume: volume must be 3D.\n" );
    }

    get_volume_sizes( volume, sizes );

    if( x < -0.5 || x > (Real) sizes[X] - 0.5 ||
        y < -0.5 || y > (Real) sizes[Y] - 0.5 ||
        z < -0.5 || z > (Real) sizes[Z] - 0.5 )
    {
        *value = get_volume_voxel_min( volume );
        if( deriv_x != (Real *) NULL )
        {
            *deriv_x = 0.0;
            *deriv_y = 0.0;
            *deriv_z = 0.0;
        }
        if( deriv_xx != (Real *) NULL )
        {
            *deriv_xx = 0.0;
            *deriv_xy = 0.0;
            *deriv_xz = 0.0;
            *deriv_yy = 0.0;
            *deriv_yz = 0.0;
            *deriv_zz = 0.0;
        }
        return;
    }

    if( x < (Real) degrees_continuity * 0.5 ||
        x > (Real) (sizes[X]-1) - (Real) degrees_continuity * 0.5 ||
        y < (Real) degrees_continuity * 0.5 ||
        y > (Real) (sizes[Y]-1) - (Real) degrees_continuity * 0.5 ||
        z < (Real) degrees_continuity * 0.5 ||
        z > (Real) (sizes[Z]-1) - (Real) degrees_continuity * 0.5 )
    {
        degrees_continuity = -1;
    }

    if( degrees_continuity < 1 && deriv_xx != (Real *) NULL )
    {
        *deriv_xx = 0.0;
        *deriv_xy = 0.0;
        *deriv_xz = 0.0;
        *deriv_yy = 0.0;
        *deriv_yz = 0.0;
        *deriv_zz = 0.0;
    }

    switch( degrees_continuity )
    {
    case -1:
        triconstant_interpolate_volume( volume, x, y, z, value,
                                        deriv_x, deriv_y, deriv_z );
        break;

    case 0:
        trilinear_interpolate_volume( volume, x, y, z,
                                      TRUE, TRUE, TRUE, value,
                                      deriv_x, deriv_y, deriv_z );
        break;

    case 1:
        triquadratic_interpolate_volume( volume, x, y, z, value,
                                         deriv_x, deriv_y, deriv_z,
                                         deriv_xx, deriv_xy, deriv_xz,
                                         deriv_yy, deriv_yz, deriv_zz );
        break;

    case 2:
        tricubic_interpolate_volume( volume, x, y, z, value,
                                     deriv_x, deriv_y, deriv_z,
                                     deriv_xx, deriv_xy, deriv_xz,
                                     deriv_yy, deriv_yz, deriv_zz );
        break;

    default:
        handle_internal_error( "evaluate_3D_volume: invalid continuity" );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_3D_slice
@INPUT      : volume
              x
              y
              z
@OUTPUT     : value
              deriv_x
              deriv_y
              deriv_xx
              deriv_xy
              deriv_yy
@RETURNS    : 
@DESCRIPTION: Takes a voxel position and evaluates the value within
              the volume by bilinear interpolation within the slice.
              If deriv_x is not a null pointer, then the
              derivatives are passed back.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

public  void   evaluate_3D_slice(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    Real           *value,
    Real           *deriv_x,
    Real           *deriv_y,
    Real           *deriv_xx,
    Real           *deriv_xy,
    Real           *deriv_yy )
{
    int           sizes[MAX_DIMENSIONS];

    if( get_volume_n_dimensions(volume) != 3 )
    {
        handle_internal_error( "evaluate_3D_slice: volume must be 3D.\n" );
    }

    get_volume_sizes( volume, sizes );

    if( x < 1.0 || x > (Real) sizes[X] - 2.0 ||
        y < 1.0 || y > (Real) sizes[Y] - 2.0 ||
        z < 1.0 || z > (Real) sizes[Z] - 2.0 )
    {
        *value = get_volume_voxel_min( volume );
        if( deriv_x != (Real *) NULL )
        {
            *deriv_x = 0.0;
            *deriv_y = 0.0;
        }
        if( deriv_xx != (Real *) NULL )
        {
            *deriv_xx = 0.0;
            *deriv_xy = 0.0;
            *deriv_yy = 0.0;
        }
        return;
    }

    bicubic_interpolate_volume( volume, x, y, z, value,
                                deriv_x, deriv_y,
                                deriv_xx, deriv_xy, deriv_yy );
}
