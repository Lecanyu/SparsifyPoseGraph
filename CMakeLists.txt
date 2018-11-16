CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(sparsifier)

#SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-multichar -std=c++0x")

IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
ENDIF()

#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")   

# Canyu Le, Close Verbose flag to refrain meaningless log
# --------------------
ADD_DEFINITIONS(-DG2S_MANIFOLD_G2O -DG2S_QUATERNIONS )
# --------------------

#ADD_DEFINITIONS(-DG2S_MANIFOLD_G2O -DG2S_QUATERNIONS -DG2S_VERBOSE)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(G2O REQUIRED)
FIND_PACKAGE(ISAM REQUIRED)
FIND_PACKAGE(Cholmod REQUIRED)
FIND_PACKAGE(Threads REQUIRED)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    FIND_PACKAGE(BLAS REQUIRED)
    SET(EXTRA_LIBRARIES
        ${BLAS_Accelerate_LIBRARY}
        ${COLAMD_LIBRARY}
        ${AMD_LIBRARY})
ENDIF()

#IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
#	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#ENDIF()
#IF((${CMAKE_BUILD_TYPE} MATCHES "Debug") AND NOT (${G2O_CORE_LIBRARY_DEBUG} MATCHES ".*NOTFOUND$"))
#    SET(G2O_LIBRARIES
#        ${G2O_CORE_LIBRARY_DEBUG}
#        ${G2O_STUFF_LIBRARY_DEBUG}
#        ${G2O_TYPES_SLAM2D_DEBUG}
#        ${G2O_TYPES_SLAM3D_DEBUG}
#        ${G2O_SOLVER_CHOLMOD_DEBUG})
#ELSE()
#    SET(G2O_LIBRARIES
#        ${G2O_CORE_LIBRARY}
#        ${G2O_STUFF_LIBRARY}
#        ${G2O_TYPES_SLAM2D}
#        ${G2O_TYPES_SLAM3D}
#        ${G2O_SOLVER_CHOLMOD})
#ENDIF()

SET(G2O_LIBRARIES
        ${G2O_CORE_LIBRARY}
        ${G2O_STUFF_LIBRARY}
        ${G2O_TYPES_SLAM2D}
        ${G2O_TYPES_SLAM3D}
        ${G2O_SOLVER_CHOLMOD})

SET(CMAKE_ECLIPSE_MAKE_ARGUMENTS "-j")

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ISAM_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${G2O_INCLUDE_DIR})

#ADD_DEFINITIONS(-DG2S_QUATERNIONS)

ADD_DEFINITIONS(-DDLIB_NO_GUI_SUPPORT)

SET(common_SOURCES
    src/pqn/pqn_optimizer.cpp
    src/pqn/line_search.cpp
    src/pqn/pqn_general.cpp
    src/isometryxd.cpp
    src/logdet_function.cpp
    src/optimizer.cpp
    src/utils.cpp
    src/type_info.cpp
    src/edge_types.cpp
    src/vertex_types.cpp
    src/glc_edge.cpp
    src/pseudo_chow_liu.cpp
    src/evaluate.cpp
    src/decimation.cpp
    src/compute_substitute_edge.cpp
    src/topology_provider_glc.cpp
    src/vertex_remover.cpp
    src/graph_wrapper_g2o.cpp
    src/graph_wrapper_isam.cpp)
    
ADD_LIBRARY(common STATIC ${common_SOURCES})

# Declare executables
ADD_EXECUTABLE(sparsifier src/main.cpp)
ADD_EXECUTABLE(kld_compare src/kld_compare.cpp)
ADD_EXECUTABLE(test_marginalize_se3 src/test_marginalize_se3.cpp)
ADD_EXECUTABLE(test_multi_edge src/test_multi_edge.cpp)
ADD_EXECUTABLE(test_glc_kld src/test_glc_kld.cpp)
ADD_EXECUTABLE(test_logdet src/test_logdet.cpp)
ADD_EXECUTABLE(test_marginalize_within_window_2 src/test_marginalize_within_window.cpp)

TARGET_LINK_LIBRARIES(common
    ${CMAKE_THREAD_LIBS_INIT}
    ${ISAM_LIBRARY}
    ${EXTRA_LIBRARIES}
    ${CHOLMOD_LIBRARY}
    ${G2O_LIBRARIES})
    
TARGET_LINK_LIBRARIES(sparsifier common)
TARGET_LINK_LIBRARIES(kld_compare common)
TARGET_LINK_LIBRARIES(test_marginalize_se3 common)
TARGET_LINK_LIBRARIES(test_multi_edge common)
TARGET_LINK_LIBRARIES(test_glc_kld common)
TARGET_LINK_LIBRARIES(test_logdet common)
TARGET_LINK_LIBRARIES(test_marginalize_within_window_2 common)

# Canyu Le
# Add boost library
# --------------------------
set(BOOST_LIB_DIR /usr/lib/x86_64-linux-gnu)
file(GLOB BOOST_LIBS ${BOOST_LIB_DIR}/libboost*)
target_link_libraries(test_marginalize_within_window_2 ${BOOST_LIBS})
# --------------------------