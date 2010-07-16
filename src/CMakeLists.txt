CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
INCLUDE_DIRECTORIES(. ../inc)

FIND_LIBRARY(JPEG NAMES jpeg)

SET(files 
    camellia.cpp 
    cam_3d.c 
    cam_project_3d_to_2d.c
    cam_3d_viewer.c
    cam_arithmetics.c 
    cam_capture.cpp 
    cam_draw.c 
    cam_keypoints.c 
    cam_keypoints_alt.c
    cam_keypoints_recurs.c 
    cam_keypoints_descriptor.c 
    cam_keypoints_matching.c 
    cam_keypoints_tracker.c 
    cam_keypoints_tracker2.c 
    cam_harris.c 
    cam_histogram.c 
    cam_hls.c 
    cam_hough.c 
    cam_hvsumming.c 
    cam_integralimage.c 
    cam_io.c 
    cam_linear_filtering.c 
    cam_LUT.c 
    cam_ME.c 
    cam_measures.c 
    cam_median_filtering.c 
    cam_morphomaths.c 
    cam_RLE_labelling.c 
    cam_RLE_morpho.c 
    cam_RLE_utils.c 
    cam_SAD.c 
    cam_utils.c 
    cam_warping.c 
    cam_watershed.c 
    cam_yuv.c 
    cam_error.c
)

ADD_LIBRARY(camellia ${files})
IF(JPEG)
    ADD_DEFINITIONS(-DHAVE_LIBJPEG)
    TARGET_LINK_LIBRARIES(camellia ${JPEG})
ENDIF(JPEG)
#ADD_LIBRARY(camellia SHARED ${files})

INSTALL(TARGETS camellia ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
INSTALL(FILES ../inc/camellia.h DESTINATION include)
