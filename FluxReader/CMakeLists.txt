# cmake .
#
# Required environment variables:
#   ROOTSYS  e.g. /nusoft/app/externals/root/v5_34_18a/Linux64bit+2.6-2.5-nu-e5-debug
#   DK2NU    e.g. /nusoft/app/externals/dk2nu/v01_00_00/Linux64bit+2.6-2.5/debug-e5/dk2nu
#
#------------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FluxReader)

set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                      $ENV{ROOTSYS}/cmake/modules
                      $ENV{ROOTSYS}/../source/root/cmake/modules
                      $ENV{ROOTSYS}/etc/cmake 
                      ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_FLAGS "-std=c++0x")

find_package(ROOT REQUIRED)

set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                        ${PROJECT_SOURCE_DIR}/include
                        $ENV{ROOTSYS}/include
                        $ENV{DK2NU}/tree)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${LINK_DIRECTORIES}
                     ${CMAKE_CURRENT_BINARY_DIR}/lib
                     $ENV{DK2NU}/lib)

link_directories(${LINK_DIRECTORIES})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project

include_directories(${PROJECT_SOURCE_DIR}/include
                    $ENV{ROOTSYS}/include
                    $ENV{DK2NU}/tree)

set(_root_internal_libraries
    Core
    Cint
    RIO
    Net
    Hist
    Graf
    Graf3d
    Gpad
    Tree
    Rint
    Postscript
    Matrix
    Physics
    MathCore
    Thread
    EG
)
set(Root_LIBRARIES ${_root_internal_libraries})

link_directories($ENV{ROOTSYS}/lib)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB FluxReader_SRCS ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB FluxReader_HDRS ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(FluxReader SHARED ${FluxReader_SRCS})
target_link_libraries(FluxReader ${ROOT_LIBRARIES} -lPhysics dk2nuTree)

install(TARGETS FluxReader DESTINATION lib)


