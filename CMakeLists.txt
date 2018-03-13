cmake_minimum_required(VERSION 3.1)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-pragmas")


# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   OFF)

# libigl options: choose your dependencies (by default everything is OFF except opengl)
if (USE_OSMESA)
    option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  OFF)
    option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           OFF)
    option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"    OFF)
else ()
    option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
    option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
    option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"    ON)
endif ()

option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

# Adding libigl: choose the path to your local copy libigl
# This is going to compile everything you requested
find_package(LIBIGL REQUIRED)


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

find_package(HDF5 REQUIRED)
if (NOT HDF5_FOUND)
  message(FATAL_ERROR "HDF5 not found")
endif()
include_directories(${HDF5_INCLUDE_DIRS})

find_package(Threads REQUIRED)

# Add project files
add_executable(${PROJECT_NAME} ${SRC_FILES} ${RENDERER_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::png)

if (USE_OSMESA)
    target_link_libraries(${PROJECT_NAME} ${OSMESA_LIBRARY})
else()
    target_link_libraries(${PROJECT_NAME} igl::opengl_glfw)
endif()

target_link_libraries(${PROJECT_NAME} ${GFLAGS_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${GLOG_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES})

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(${PROJECT_NAME} PUBLIC "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
  target_link_libraries(${PROJECT_NAME} "${CMAKE_THREAD_LIBS_INIT}")
endif()
