#TODO convert testing routines from Makefile.am here

STRING(REPLACE ";" " " MY_ENV "${MINC_TEST_ENVIRONMENT}")

set_directory_properties( PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")

configure_file(minc_wrapper.cmake ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper @ONLY)

MACRO(add_minc_test name command )
 ADD_TEST( ${name}
    cmake -E chdir ${CMAKE_CURRENT_BINARY_DIR} ${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(add_minc_test)

SET(ENV{DYLD_LIBRARY_PATH} "/Users/vfonov/src/build/minc-toolkit-itk4/external/opt/minc-itk4/lib")

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc
  DEPENDS mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper make_phantom -clobber -ellipse -nele 64 64 64 -step 2 2 2 -start -64 -64 -64 -center 0 0 0 ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/test1.xfm
  DEPENDS mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper param2xfm -clobber -translation 10 20 30 ${CMAKE_CURRENT_BINARY_DIR}/test1.xfm
)
# This could be done by applying test1.xfm to ellipse0,
# but we choose to do it by making a phantom with a different centre.
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc
  DEPENDS mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper make_phantom -clobber -ellipse -nele 64 64 64 -step 2 2 2 -start -64 -64 -64 -center 10 20 30 ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc
)

add_custom_command(
  OUTPUT   ${CMAKE_CURRENT_BINARY_DIR}/test2.xfm
  COMMAND  ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10 ${CMAKE_CURRENT_BINARY_DIR}/test2.xfm
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse2.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/test2.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_BINARY_DIR}/test2.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse2.mnc
)
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/test3.xfm
  DEPENDS mni_autoreg  
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10  -scales 1.7 1.7 1.7 ${CMAKE_CURRENT_BINARY_DIR}/test3.xfm
)
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse3.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/test3.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_BINARY_DIR}/test3.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc  ${CMAKE_CURRENT_BINARY_DIR}/ellipse3.mnc
)
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincmath ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp1.mnc -mult -constant 4.5  -clob
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_SOURCE_DIR}/tps.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp1.mnc ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4.mnc 
  COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp1.mnc
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse1_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse1
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse2_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse2.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/ellipse2.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse2
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse3_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse3.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/ellipse3.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse3
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4.mnc ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse4.mnc 
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -xdircos 0.891006529493696 0 0.453990489327253 -ydircos 0 1 0 -zdircos -0.453990489327253 0 0.891006529493696 ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse4.mnc 
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse4_dxyz.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse4.mnc mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/ellipse4.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse4
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm
  DEPENDS mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper param2xfm -clobber -translation 3 2 0 ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_x.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp2.mnc
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincreshape -clobber -dimrange xspace=32,1 ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp2.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_x.mnc
  COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp2.mnc
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_y.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp3.mnc
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincreshape -clobber -dimrange yspace=35,1 ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp3.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_y.mnc
  COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp3.mnc
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_z.mnc
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample -clobber -transformation ${CMAKE_CURRENT_BINARY_DIR}/test4.xfm -like ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp4.mnc
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincreshape -clobber -dimrange zspace=35,1 ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp4.mnc ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_z.mnc
  COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/ellipse_tmp4.mnc
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc 
  DEPENDS mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper make_phantom -nele 64 64 64 -step 2 2 2 -start -64 -64 -64 -ellipse -center 0 0 0 -clobber ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc
)

add_custom_command(
 OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/test.xfm 
 DEPENDS mni_autoreg
 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper param2xfm -center 0 0 0 -rot -4 7 10 -trans 5 2 -6 -clobber ${CMAKE_CURRENT_BINARY_DIR}/test.xfm
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/object2.mnc 
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc  ${CMAKE_CURRENT_BINARY_DIR}/test.xfm mni_autoreg
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincresample ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc ${CMAKE_CURRENT_BINARY_DIR}/object2.mnc -like ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc -transformation ${CMAKE_CURRENT_BINARY_DIR}/test.xfm -clobber
)

add_custom_command(
 OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/object1_dxyz.mnc 
 DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc mni_autoreg
 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc ${CMAKE_CURRENT_BINARY_DIR}/object1
)
 
add_custom_command(
 OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/object2_dxyz.mnc
 DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/object2.mnc mni_autoreg
 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/minc_wrapper mincblur -clobber -gradient -fwhm 6 ${CMAKE_CURRENT_BINARY_DIR}/object2.mnc ${CMAKE_CURRENT_BINARY_DIR}/object2
)

add_custom_target(create_test_data ALL 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse1.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse2.mnc
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse3.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse4.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse1_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse2_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse3_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse4_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/warped_ellipse4_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_x.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_y.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/ellipse0_slice_z.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/test1.xfm 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/test2.xfm 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/test3.xfm 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/object1.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/object2.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/object1_dxyz.mnc 
 DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/object2_dxyz.mnc
 DEPENDS mni_autoreg
)

add_minc_test(param2xfm           ${CMAKE_CURRENT_SOURCE_DIR}/param2xfm.test.cmake)
add_minc_test(minctracc_linear    ${CMAKE_CURRENT_SOURCE_DIR}/minctracc.test1.cmake)
add_minc_test(minctracc_nonlinear ${CMAKE_CURRENT_SOURCE_DIR}/minctracc.test2.cmake)
