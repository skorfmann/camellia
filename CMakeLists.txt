CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Camellia)

SUBDIRS(src)

IF(CMAKE_HOST_WIN32)
    SET(libs winmm)    
ENDIF(CMAKE_HOST_WIN32)

INCLUDE_DIRECTORIES(inc)
ADD_EXECUTABLE(cam_demo cam_demo.c cam_demo_cpp.cpp)
TARGET_LINK_LIBRARIES(cam_demo camellia ${libs})

ADD_EXECUTABLE(cam_myobjects cam_myobjects.c)
TARGET_LINK_LIBRARIES(cam_myobjects camellia ${libs} m)

