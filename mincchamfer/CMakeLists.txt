ADD_EXECUTABLE(mincchamfer	
 chamfer.c chamfer.h  
 mincchamfer.c mincchamfer.h)

TARGET_LINK_LIBRARIES(mincchamfer Proglib)


INSTALL(TARGETS 
  mincchamfer 
  DESTINATION bin)
