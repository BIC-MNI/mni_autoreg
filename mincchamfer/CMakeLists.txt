ADD_EXECUTABLE(mincchamfer	
 chamfer.c chamfer.h  
	mincchamfer.c mincchamfer.h)

TARGET_LINK_LIBRARIES(mincchamfer Proglib)