#TODO convert testing routines from Makefile.am here

MACRO(add_minc_test name command )
 ADD_TEST( ${name}
    ${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(add_minc_test)


add_custom_command(
  OUTPUT  ellipse0.mnc
  COMMAND make_phantom -clobber -ellipse -nele 64 64 64 -step 2 2 2 -start -64 -64 -64 -center 0 0 0 ellipse0.mnc 
)

add_custom_command(
  OUTPUT test1.xfm
  COMMAND param2xfm -clobber -translation 10 20 30 test1.xfm
)
# This could be done by applying test1.xfm to ellipse0,
# but we choose to do it by making a phantom with a different centre.
add_custom_command(
  OUTPUT ellipse1.mnc
  COMMAND make_phantom -clobber -ellipse -nele 64 64 64 -step 2 2 2 -start -64 -64 -64 -center 10 20 30 ellipse1.mnc
)

add_custom_command(
  OUTPUT test2.xfm
  COMMAND  param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10 test2.xfm
)

add_custom_command(
  OUTPUT ellipse2.mnc
  DEPENDS ellipse0.mnc test2.xfm
  COMMAND mincresample -clobber -transformation test2.xfm -like ellipse0.mnc ellipse0.mnc ellipse2.mnc
)
add_custom_command(
  OUTPUT test3.xfm
  COMMAND param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10  -scales 1.7 1.7 1.7 test3.xfm
)
add_custom_command(
  OUTPUT ellipse3.mnc
  DEPENDS ellipse0.mnc test3.xfm
  COMMAND mincresample -clobber -transformation test3.xfm -like ellipse0.mnc ellipse0.mnc  ellipse3.mnc
)
add_custom_command(
  OUTPUT warped_ellipse4.mnc
  DEPENDS ellipse0.mnc
  COMMAND mincmath ellipse0.mnc ellipse_tmp1.mnc -mult -constant 4.5  -clob
  COMMAND mincresample -clobber -transformation ${CMAKE_CURRENT_SOURCE_DIR}/tps.xfm -like ellipse0.mnc ellipse_tmp1.mnc warped_ellipse4.mnc 
  COMMAND rm -f ellipse_tmp1.mnc
)

add_custom_command(
  OUTPUT ellipse0_dxyz.mnc
  DEPENDS ellipse0.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 ellipse0.mnc ellipse0
)

add_custom_command(
  OUTPUT ellipse1_dxyz.mnc
  DEPENDS ellipse1.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 ellipse1.mnc ellipse1
)

add_custom_command(
  OUTPUT ellipse2_dxyz.mnc
  DEPENDS ellipse2.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 ellipse2.mnc ellipse2
)

add_custom_command(
  OUTPUT ellipse3_dxyz.mnc
  DEPENDS ellipse3.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 ellipse3.mnc ellipse3
)

add_custom_command(
  OUTPUT warped_ellipse4_dxyz.mnc
  DEPENDS warped_ellipse4.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 warped_ellipse4.mnc warped_ellipse4
)

add_custom_command(
  OUTPUT ellipse4.mnc 
  DEPENDS ellipse0.mnc
  COMMAND mincresample -clobber -xdircos 0.891006529493696 0 0.453990489327253 -ydircos 0 1 0 -zdircos -0.453990489327253 0 0.891006529493696 ellipse0.mnc ellipse4.mnc 
)

add_custom_command(
  OUTPUT ellipse4_dxyz.mnc
  DEPENDS ellipse4.mnc
  COMMAND mincblur -clobber -gradient -fwhm 6 ellipse4.mnc ellipse4
)

add_custom_command(
  OUTPUT test4.xfm
  COMMAND param2xfm -clobber -translation 3 2 0 test4.xfm
)

add_custom_command(
  OUTPUT ellipse0_slice_x.mnc
  DEPENDS ellipse0.mnc test4.xfm
  COMMAND mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp2.mnc
  COMMAND mincreshape -clobber -dimrange xspace=32,1 ellipse_tmp2.mnc ellipse0_slice_x.mnc
  COMMAND rm -f ellipse_tmp2.mnc
)

add_custom_command(
  OUTPUT ellipse0_slice_y.mnc
  DEPENDS ellipse1.mnc test4.xfm
  COMMAND mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp3.mnc
  COMMAND mincreshape -clobber -dimrange yspace=35,1 ellipse_tmp3.mnc ellipse0_slice_y.mnc
  COMMAND rm -f ellipse_tmp3.mnc
)

add_custom_command(
  OUTPUT ellipse0_slice_z.mnc
  DEPENDS ellipse1.mnc test4.xfm
  COMMAND mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp4.mnc
  COMMAND mincreshape -clobber -dimrange zspace=35,1 ellipse_tmp4.mnc ellipse0_slice_z.mnc
  COMMAND rm -f ellipse_tmp4.mnc
)

add_custom_target(create_test_data ALL echo 
 DEPENDS  ellipse0.mnc ellipse1.mnc ellipse2.mnc ellipse3.mnc 
 DEPENDS  ellipse4.mnc warped_ellipse4.mnc ellipse0_dxyz.mnc 
 DEPENDS  ellipse1_dxyz.mnc ellipse2_dxyz.mnc ellipse3_dxyz.mnc 
 DEPENDS  ellipse4_dxyz.mnc warped_ellipse4_dxyz.mnc ellipse0_slice_x.mnc 
 DEPENDS  ellipse0_slice_y.mnc ellipse0_slice_z.mnc test1.xfm 
 DEPENDS  test2.xfm test3.xfm )

