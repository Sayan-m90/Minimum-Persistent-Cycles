#project( pers2cyc_fin )

cmake_minimum_required( VERSION 3.5 FATAL_ERROR )
#cmake_policy( SET CMP0017 NEW )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3")


include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" )
link_directories(${CMAKE_SOURCE_DIR}/src/lib)
# Boost
find_package( Boost 1.5 COMPONENTS filesystem program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

#GUDHI
set(GUDHI_INCLUDE_DIR /usr/local/include)
include_directories(${GUDHI_INCLUDE_DIR})

include_directories( "." ".." ) 

add_executable( pers2cyc_fin
	twocyc.cpp
	Legal.h
	ParseCommand.h
	Graph.h
	) 

target_link_libraries( pers2cyc_fin ${Boost_LIBRARIES})
