# --------------------------------------------------------------------
# 
# Makefile to build tests
#



clean:
	rm -f *.mnc *.xfm *.obj

test.xfm: 
	param2xfm -center 0 0 0 -rot -4 7 10 -trans 5 2 -6 \
		-clobber test.xfm 

object1.mnc: 
	make_phantom -nele 64 64 64 -step 2 2 2 \
		-start -64 -64 -64 -ellipse -center 0 0 0 \
		-clobber object1.mnc

object2.mnc: object1.mnc test.xfm
	mincresample object1.mnc object2.mnc -like object1.mnc -transformation test.xfm -clobber

            
object1_dxyz.mnc: object1.mnc
	mincblur -clobber -gradient -fwhm 6 object1.mnc object1

object2_dxyz.mnc: object2.mnc
	mincblur -clobber -gradient -fwhm 6 object2.mnc object2

object1_blur.mnc: object1.mnc
	mincblur -clobber -gradient -fwhm 6 object1.mnc object1

object2_blur.mnc: object2.mnc
	mincblur -clobber -gradient -fwhm 6 object2.mnc object2

o1.mnc: object1_blur.mnc
	autocrop object1_blur.mnc -isostep 3 o1.mnc

o2.mnc: object2_blur.mnc
	autocrop object2_blur.mnc -isostep 3 o2.mnc

minctracc.test: object1_dxyz.mnc object2_dxyz.mnc
	../animal -identity object1_dxyz.mnc object2_dxyz.mnc \
		  -est_center -debug -simplex 10 -lsq6 -step 8 8 8 \
		 -clobber output.xfm 
	@echo ""
	@echo "Minctracc should have found rotations close to -4,7,10;"
	@echo "and translations  close to 5,2,-6"
	xfm2param -center 0 0 0 output.xfm
	touch $@
	@echo "minctracc successful"
	@echo ""



#
# test animal with 5 iterations, simplex optimization
# and a small offset between objects

def5: object1_dxyz.mnc object2_dxyz.mnc
	param2xfm -clob -translation 5 2 -6 trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_dxyz.mnc object2_dxyz.mnc def1.xfm \
	-est_center -step 6 6 6 -tol 0.000001 \
	-transformation trans.xfm \
	-iterations 5 -weight 0.5 -stiff 0.5 -similar 0.5

trans.xfm:
	param2xfm -clob -translation 5 2 -6 trans.xfm

defq: object1_dxyz.mnc object2_dxyz.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_dxyz.mnc object2_dxyz.mnc def_q1.xfm \
	-est_center -step 6 6 6 -tol 0.000001 \
	-transformation trans.xfm -use_local -quad \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5

defqs: object1_blur.mnc object1_blur.mnc
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_qs.xfm \
	-est_center -step 6 6 2 -tol 0.000001 \
	-identity -use_noniso -use_local -quad \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5

defq5: object1_dxyz.mnc object2_dxyz.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_dxyz.mnc object2_dxyz.mnc def_q5.xfm \
	-est_center -step 6 6 6 -tol 0.000001 \
	-transformation trans.xfm -use_local -quad \
	-iterations 5 -weight 1.0 -stiff 0.5 -similar 0.5 

mask.mnc: object1.mnc object2.mnc
	mincmath -add object{1,2}.mnc mask.mnc

defqn: mask.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_qn1.xfm \
	-est_center -step 6 6 6 -tol 0.000001 \
	-transformation trans.xfm -use_local -quad \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso
	def_to_lines pr_axes2_0.xfm pr_qn1.obj mask.mnc .3

#
# map a simple object to itself to ensure that the eigen values returned
# point in the correct direction.
#

defqn1: object1.mnc object1_blur.mnc 
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object1_blur.mnc def_qn1.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -quad \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso
	def_to_lines pr_axes2_0.xfm pr_qn1.obj object1.mnc .3

defqn4: object1.mnc object1_blur.mnc 
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object1_blur.mnc def_qn4.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -quad \
	-iterations 4 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso
	def_to_lines pr_axes2_0.xfm pr_qn4.obj object1.mnc .3

#
#  now try to map one object onto another.  (max disp is 1cm)
#
#
defqn4_2: mask.mnc object1_blur.mnc object2_blur.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_qn4_2.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -quad -transformation trans.xfm \
	-iterations 4 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso
	def_to_lines pr_axes2_0.xfm pr_qn4_2.obj mask.mnc .3

defqn1_2: mask.mnc object1_blur.mnc object2_blur.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_qn1_2.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -quad -transformation trans.xfm \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso

defqn2_2: mask.mnc object1_blur.mnc object2_blur.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_qn2_2.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -quad -transformation trans.xfm \
	-iterations 2 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso

defs1_2: mask.mnc object1_blur.mnc object2_blur.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_blur.mnc object2_blur.mnc def_s1_2.xfm \
	-est_center -step 4 4 4 -tol 0.000001 \
	 -use_local -transformation trans.xfm \
	-iterations 1 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso


defn5: object1_dxyz.mnc object2_dxyz.mnc trans.xfm
	../animal -debug -xcorr -nonlin -clobber -verbose 2 \
	object1_dxyz.mnc object2_dxyz.mnc def_qn5.xfm \
	-est_center -step 6 6 6 -tol 0.000001 \
	-transformation trans.xfm -use_local -quad \
	-iterations 5 -weight 1.0 -stiff 0.5 -similar 0.5 \
	-use_noniso


