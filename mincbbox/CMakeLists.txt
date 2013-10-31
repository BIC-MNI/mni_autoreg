ADD_EXECUTABLE(mincbbox  mincbbox.c mincbbox.h)
TARGET_LINK_LIBRARIES(mincbbox Proglib)

INSTALL(TARGETS 
  mincbbox 
  DESTINATION bin)
