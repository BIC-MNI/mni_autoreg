ADD_EXECUTABLE(mincblur
	apodize_data.c 
	blur_support.c blur_support.h 
	blur_volume.c blur_volume.h 
	fft.c 
	gradient_volume.c 
	gradmag_volume.c gradmag_volume.h 
	kernel.h 
	mincblur.c mincblur.h)

TARGET_LINK_LIBRARIES(mincblur Proglib)


INSTALL(TARGETS 
  mincblur 
  DESTINATION bin)