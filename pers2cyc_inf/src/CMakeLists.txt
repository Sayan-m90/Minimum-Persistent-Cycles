cmake_minimum_required( VERSION 3.5 FATAL_ERROR )
#cmake_policy( SET CMP0017 NEW )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3" )

message( "Compiler: ${CMAKE_CXX_COMPILER}" )
message( "Option: ${CMAKE_CXX_FLAGS}" )

# preprocessor macros
#add_compile_definitions( DEBUG_OPT )

include_directories( "." )

# Boost
find_package( Boost 1.5 COMPONENTS filesystem program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} )

# CGAL
find_package( CGAL )

# Eigen
include_directories( $ENV{Eigen_INC_DIR} )

# GUDHI
include_directories( $ENV{GUDHI_INC_DIR} )

add_executable( pers2cyc_inf
    main.cpp
    cell_complex.cpp
    cube_complex_bm.cpp
    utils.cpp
    geom_utils.cpp
    tests.cpp
    inf_pers2cyc.cpp
    mesh_writer.cpp
	)

target_link_libraries( pers2cyc_inf 
    ${Boost_PROGRAM_OPTIONS_LIBRARY} 
    CGAL::CGAL )
