cmake_minimum_required(VERSION 2.6)

# Set whether use OSMesa
set(USE_OSMESA CACHE STRING FALSE)

if (USE_OSMESA)
    project(OSMesaRenderer)
    message(STATUS "== Use OSMesa Mesh Renderer ==")
else()
    project(LibiglViewer)
    message(STATUS "== Use Libigl Mesh Renderer ==")
endif()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(LIBIGL)

if (NOT LIBIGL_FOUND)
    message(FATAL_ERROR "libigl not found --- You can download it using: \n git clone --recursive https://github.com/libigl/libigl.git ${PROJECT_SOURCE_DIR}/../libigl")
endif()

# Compilation flags: adapt to your needs
if(MSVC)
    # Enable parallel compilation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
    # Libigl requires a modern C++ compiler that supports c++11
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   OFF)

# libigl options: choose your dependencies (by default everything is OFF except opengl)
if (USE_OSMESA)
    option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  OFF)
    option(LIBIGL_WITH_GLFW        "Use GLFW"           OFF)
else ()
    option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
    option(LIBIGL_WITH_GLFW        "Use GLFW"           ON)
endif ()

option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
    find_package(CGAL REQUIRED)
    set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
    include(${CGAL_USE_FILE})
endif()

# Adding libigl: choose the path to your local copy libigl
# This is going to compile everything you requested
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../shared/cmake" "libigl")

# libigl information
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

# Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

include_directories(
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/base"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")

FILE(GLOB SRC_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/base/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/*/*.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/*/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/utils/*/*.cpp")

if (USE_OSMESA)
    add_definitions(-DUSE_OSMESA)
    find_package(OSMesa REQUIRED)
    if (NOT OSMESA_FOUND)
        message(FATAL_ERROR "OSMesa not found")
    endif()
    include_directories(${OSMESA_INCLUDE_DIR})
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/renderers/osmesa")
    FILE(GLOB RENDERER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/renderers/osmesa/*.cc")
else()
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/renderers/libigl")
    FILE(GLOB RENDERER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/renderers/libigl/*.cc")
endif()

# Add libraries
find_package(GFlags REQUIRED)
if (NOT GFLAGS_FOUND)
    message(FATAL_ERROR "gflags not found")
endif()
include_directories(${GFLAGS_INCLUDE_DIRS})

find_package(Glog REQUIRED)
if (NOT GLOG_FOUND)
    message(FATAL_ERROR "glog not found")
endif()
include_directories(${GLOG_INCLUDE_DIRS})

# Add project files
add_executable(${PROJECT_NAME} ${SRC_FILES} ${RENDERER_FILES} ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

if (USE_OSMESA)
    target_link_libraries(${PROJECT_NAME} ${OSMESA_LIBRARY})
endif()

target_link_libraries(${PROJECT_NAME} ${GFLAGS_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${GLOG_LIBRARIES})
