cmake_minimum_required(VERSION 2.6)
project (mih)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "DEBUG")
ENDIF()

find_library (HDF5_LIB hdf5)


set (CMAKE_CXX_FLAGS "-DSINGLE_CORE -O2 -fopenmp -DNDEBUG")
# -Wall -march=native -mtune=native")

include_directories ("${PROJECT_SOURCE_DIR}/include")

add_executable (mih src/array32.cpp src/mihasher.cpp src/sparse_hashtable.cpp src/bucket_group.cpp interface/mih_interface.cpp src/reorder.cpp interface/loadVar.cpp interface/saveRes.cpp)
target_link_libraries (mih ${HDF5_LIB})

project (linscan)
add_executable (linscan interface/linscan_interface.cpp src/linscan.cpp interface/loadVar.cpp interface/saveRes.cpp)
target_link_libraries (linscan ${HDF5_LIB})

FIND_PACKAGE(Boost 1.45.0)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" "${PROJECT_SOURCE_DIR}/include" "/usr/include/python2.7")
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

  ADD_LIBRARY(mihmodule SHARED src/mihasher.cpp src/mihmodule.cpp src/array32.cpp src/sparse_hashtable.cpp src/bucket_group.cpp interface/mih_interface.cpp src/reorder.cpp interface/loadVar.cpp interface/saveRes.cpp)
  TARGET_LINK_LIBRARIES(mihmodule ${HDF5_LIB} ${Boost_LIBRARIES} ) 
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()
