cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

SET( CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-long-long")

#---------------------------
set( PackageName CaloGeo )
#---------------------------

project(${PackageName})

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 0 )
SET( ${PackageName}_VERSION_PATCH 0 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${CMAKE_SOURCE_DIR}/cmake ) 
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH  
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..." 
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



#---Options-------------------------------------------------------------------------

option(BUILD_TESTING "Enable and build tests" ON)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4)
find_package(ROOT REQUIRED)
#find_package(Geant4 REQUIRED ) 
find_package( Boost REQUIRED ) 

set( ROOT_COMPONENT_LIBRARIES Geom Reflex)

include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
include_directories( ${CMAKE_SOURCE_DIR}/detector/include )
include_directories( 
  SYSTEM
  ${DD4hep_INCLUDE_DIRS}
  ${DD4hep_COMPONENTS_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
)

#${Geant4_INCLUDE_DIRS}

file(
  GLOB sources 
  ./SDHCAL/src/*.cpp 
)

#Create target, library, rootmap, install
add_dd4hep_plugin(${PackageName} SHARED ${sources})
target_link_libraries(${PackageName} ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} 
# ${Geant4_LIBRARIES}
)

#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})

