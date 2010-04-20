CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Camellia)

SUBDIRS(src)

INCLUDE_DIRECTORIES(inc)
ADD_EXECUTABLE(cam_demo cam_demo.c cam_demo_cpp.cpp)
TARGET_LINK_LIBRARIES(cam_demo camellia)

