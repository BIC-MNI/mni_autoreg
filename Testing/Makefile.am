# --------------------------------------------------------------------
# 
# Makefile to build tests
#
TEST_LIST = \
	param2xfm.test \
	xfm2param.test \
	make_phantom.test \
	mincresample.test \
	mincblur.test \
	xcorr_vol.test \
	volume_cog.test \
	minctracc.test



check-local: banner $(TEST_LIST)
	@echo "Tests completed successfully."

clean:
	rm -f *.test *.act *.mnc *.xfm *reals

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

