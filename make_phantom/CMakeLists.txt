ADD_EXECUTABLE(make_phantom make_phantom.c make_phantom.h)

TARGET_LINK_LIBRARIES(make_phantom Proglib)

INSTALL(TARGETS 
  make_phantom 
  DESTINATION bin)
