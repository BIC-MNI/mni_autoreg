# Run some checks.  

# The checks use programs built elsewhere.  One must run
# 'make check' from the top-level directory, to ensure these are built.
minctraccdir = $(top_builddir)/minctracc/Main
extradir = $(top_builddir)/minctracc/Extra_progs
built_PATH = $(minctraccdir):$(extradir)

# The tests are shell scripts, so the environment ends with the
# path to the shell interpreter.
TESTS_ENVIRONMENT = PATH=$(built_PATH):$(PATH) $(SHELL)

TESTS = linear-1 linear-2 linear-3 nonlinear-2 nonlinear-3 nonlinear-4 nonlinear-5 nonlinear-6 nonlinear-7 nonlinear-8

EXTRA_DIST = $(TESTS) tps.xfm tps.tag


# Objects used in testing
aux_testfiles = \
	ellipse0.mnc \
	ellipse1.mnc \
	ellipse2.mnc \
	ellipse3.mnc \
	ellipse4.mnc \
	warped_ellipse4.mnc \
	ellipse0_dxyz.mnc \
	ellipse1_dxyz.mnc \
	ellipse2_dxyz.mnc \
	ellipse3_dxyz.mnc \
	ellipse4_dxyz.mnc \
	warped_ellipse4_dxyz.mnc \
	ellipse0_slice_x.mnc \
	ellipse0_slice_y.mnc \
	ellipse0_slice_z.mnc \
	test1.xfm \
	test2.xfm \
	test3.xfm 

check_DATA = $(aux_testfiles)

CLEANFILES = $(aux_testfiles) \
	linear-1.log linear-2.log linear-3.log \
	nonlinear-2.log nonlinear-3.log nonlinear-4.log nonlinear-5.log nonlinear-6.log nonlinear-7.log nonlinear-8.log

ellipse0.mnc: Makefile.am
	../make_phantom/make_phantom -clobber -ellipse \
	-nele 64 64 64 -step 2 2 2 -start -64 -64 -64 \
	-center 0 0 0 $@

test1.xfm: Makefile.am
	$(extradir)/param2xfm -clobber -translation 10 20 30 $@

# This could be done by applying test1.xfm to ellipse0,
# but we choose to do it by making a phantom with a different centre.
ellipse1.mnc: Makefile.am
	../make_phantom/make_phantom -clobber -ellipse \
	-nele 64 64 64 -step 2 2 2 -start -64 -64 -64 \
	-center 10 20 30 $@

test2.xfm: Makefile.am
	$(extradir)/param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10 $@

ellipse2.mnc: Makefile.am ellipse0.mnc test2.xfm
	mincresample -clobber -transformation test2.xfm \
	-like ellipse0.mnc ellipse0.mnc $@

test3.xfm: Makefile.am
	$(extradir)/param2xfm -clobber -translation 5 2 -6 -rotation -4 7 10 \
	-scales 1.7 1.7 1.7 $@

ellipse3.mnc: Makefile.am ellipse0.mnc test3.xfm
	mincresample -clobber -transformation test3.xfm \
	-like ellipse0.mnc ellipse0.mnc $@

warped_ellipse4.mnc: Makefile.am ellipse0.mnc
	mincmath ellipse0.mnc ellipse_tmp$$.mnc -mult -constant 4.5 ; \
	mincresample -clobber -transformation $(srcdir)/tps.xfm \
	-like ellipse0.mnc ellipse_tmp$$.mnc $@ ; \
	rm -f ellipse_tmp$$.mnc

ellipse0_dxyz.mnc: ellipse0.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 ellipse0.mnc ellipse0

ellipse1_dxyz.mnc: ellipse1.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 ellipse1.mnc ellipse1

ellipse2_dxyz.mnc: ellipse2.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 ellipse2.mnc ellipse2

ellipse3_dxyz.mnc: ellipse3.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 ellipse3.mnc ellipse3

warped_ellipse4_dxyz.mnc: warped_ellipse4.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 warped_ellipse4.mnc warped_ellipse4

ellipse4.mnc: Makefile.am ellipse0.mnc
	mincresample -clobber \
        -xdircos 0.891006529493696 0 0.453990489327253 \
        -ydircos 0 1 0 \
        -zdircos -0.453990489327253 0 0.891006529493696 \
	ellipse0.mnc $@

ellipse4_dxyz.mnc: ellipse4.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 ellipse4.mnc ellipse4

test4.xfm: Makefile.am
	$(extradir)/param2xfm -clobber -translation 3 2 0 $@

ellipse0_slice_x.mnc: Makefile.am ellipse0.mnc test4.xfm
	mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp$$.mnc; \
	mincreshape -clobber -dimrange xspace=32,1 ellipse_tmp$$.mnc $@

ellipse0_slice_y.mnc: Makefile.am ellipse1.mnc test4.xfm
	mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp$$.mnc; \
	mincreshape -clobber -dimrange yspace=35,1 ellipse_tmp$$.mnc $@

ellipse0_slice_z.mnc: Makefile.am ellipse1.mnc test4.xfm
	mincresample -clobber -transformation test4.xfm -like ellipse0.mnc ellipse0.mnc ellipse_tmp$$.mnc; \
	mincreshape -clobber -dimrange zspace=35,1 ellipse_tmp$$.mnc $@


# --------------------------------------------------------------------
# 
# (old) Makefile to build tests
#
TEST_LIST = \
	param2xfm.test \
	xfm2param.test \
	make_phantom.test \
	mincresample.test \
	mincblur.test \
	xcorr_vol.test \
	volume_cog.test \
	minctracc.test \
	animal.test



old-checks: banner $(TEST_LIST)
	@echo "Tests completed successfully."

clean-local:
	rm -f *.test *.act *.mnc test*.xfm output.xfm

# -------banner message------------------
banner:
	@echo "---------------------------------------------------"
	@echo "These are tests of the individual programs used in "
	@echo "the MNI linear registration package."
	@echo "---------------------------------------------------"
	@echo ""

# how to make individual objects ------------------------------

test.xfm: param2xfm.test

object1.mnc: make_phantom.test

object2.mnc: mincresample.test 

object1_dxyz.mnc: object1.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 object1.mnc object1

object2_dxyz.mnc: object2.mnc
	../mincblur/mincblur -clobber -gradient -fwhm 6 object2.mnc object2

# program tests -----------------------------------------------

mess.mincblur:
	@echo "---------------------------------------"
	@echo "testing volumetric gaussian blurring"
	@echo "(this may take a few minutes)"
	@echo ""

mincblur.test: mess.mincblur mincblur.act

mincblur.act: object1_dxyz.mnc object2_dxyz.mnc
	touch $@
	@echo "mincblur successful"
	@echo ""

mess.mincresample:
	@echo "---------------------------------------"
	@echo "testing volumetric resampling"
	@echo ""

mincresample.test: mess.mincresample mincresample.act

mincresample.act: object1.mnc test.xfm
	mincresample object1.mnc object2.mnc -like object1.mnc -transformation test.xfm -clobber
	touch $@
	@echo "mincresample successful"
	@echo ""



mess.make_phantom:
	@echo "---------------------------------------"
	@echo "testing program make_phantom"
	@echo "(this may take a few minutes)"
	@echo ""

make_phantom.test: mess.make_phantom make_phantom.act

make_phantom.act:
	../make_phantom/make_phantom -nele 64 64 64 -step 2 2 2 \
		-start -64 -64 -64 -ellipse -center 0 0 0 \
		-clobber object1.mnc
	touch $@
	@echo "make_phantom successful"
	@echo ""



mess.param2xfm:
	@echo "---------------------------------------"
	@echo "testing program param2xfm"
	@echo ""

param2xfm.test: mess.param2xfm param2xfm.act

param2xfm.act:
	../minctracc/Extra_progs/param2xfm -center 0 0 0 -rot -4 7 10 -trans 5 2 -6 \
		-clobber test.xfm 
	cat test.xfm
	touch $@
	@echo "param2xfm successful"
	@echo ""


mess.xfm2param:
	@echo "---------------------------------------"
	@echo "testing parameter extraction"
	@echo ""

xfm2param.test: mess.xfm2param xfm2param.act

xfm2param.act: test.xfm
	@echo ""
	@echo "The following should give rotations of -4,7,10;"
	@echo "translations of 5,2,-6"
	../minctracc/Extra_progs/xfm2param test.xfm
	touch $@
	@echo "xfm2param successful"
	@echo ""


mess.minctracc:
	@echo "---------------------------------------"
	@echo "testing minctracc: the registration engine"
	@echo "(this may take a few minutes)"
	@echo ""

minctracc.test: mess.minctracc minctracc.act

minctracc.act: object1_dxyz.mnc object2_dxyz.mnc
	../minctracc/Main/minctracc -identity object1_dxyz.mnc object2_dxyz.mnc \
		 -est_center -debug -simplex 10 -lsq6 -step 8 8 8 \
		-clobber output.xfm 
	@echo ""
	@echo "Minctracc should have found rotations close to -4,7,10;"
	@echo "and translations  close to 5,2,-6"
	../minctracc/Extra_progs/xfm2param -center 0 0 0 output.xfm
	touch $@
	@echo "minctracc successful"
	@echo ""


mess.animal:
	@echo "---------------------------------------"
	@echo "testing animal: the nonlinear registration engine"
	@echo "(this may take a few minutes)"
	@echo ""

animal.test: mess.animal animal.act

animal.act: object1_dxyz.mnc object2_dxyz.mnc
	../minctracc/Main/minctracc -identity object1_dxyz.mnc object2_dxyz.mnc \
		-est_center -debug  -step 10 10 10 -nonlin \
		-clobber def.xfm 
	@echo ""
	touch $@
	@echo "animal successful"
	@echo ""


mess.volumecog:
	@echo "---------------------------------------"
	@echo "testing volume center of gravity"
	@echo ""

volume_cog.test: mess.volumecog volume_cog.act

volume_cog.act: object1.mnc
	@echo "The object's center of gravity should be 0,0,0"
	../minctracc/Extra_progs/volume_cog object1.mnc
	touch $@
	@echo "volume_cog successful"
	@echo ""


mess.xcorr_vol:
	@echo "---------------------------------------"
	@echo "testing volume correlation measurement"
	@echo ""

xcorr_vol.test: mess.xcorr_vol xcorr_vol.act

xcorr_vol.act: object1.mnc object2.mnc
	@echo "The volume correlation of object1 and object2 should"
	@echo "be 0.723"
	@echo "running...."
	../minctracc/Extra_progs/xcorr_vol object1.mnc object2.mnc
	touch $@
	@echo "xcorr_vol successful"
	@echo ""

