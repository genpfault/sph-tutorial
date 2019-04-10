cmake_minimum_required( VERSION 3.10.0 )
project( sph-benchmark )
set( CMAKE_CXX_STANDARD 11 )

# external glfw
set( GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE )
set( GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE )
set( GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE )
add_subdirectory( external/glfw )

# make sure we have OpenMP C++ support
find_package( OpenMP REQUIRED CXX )

# external includes
include_directories( SYSTEM "external/glm" )
include_directories( SYSTEM "external/glad/include" )

# application code
include_directories( "include" )
add_executable(
    sph-benchmark
    "src/main.cpp"
    "external/glad/src/glad.c" )
target_link_libraries(
    sph-benchmark
    PRIVATE glfw
    PRIVATE OpenMP::OpenMP_CXX )
