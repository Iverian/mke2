cmake_minimum_required(VERSION 3.0)
project(mke2 VERSION 0.1.0 LANGUAGES CXX)

###################################################
# Include additional cmake packages
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
###################################################

set(CMAKE_CXX_STANDARD 17)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

###################################################
# Declare dependencies
# find_package(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
###################################################

###################################################
# Find source files
file(GLOB MAIN_HEADERS "include/*.hpp")

set(MAIN_EXECUTABLE "src/main.cpp")
set(MAIN_SOURCES "")
file(GLOB _SRC_CHILDREN "src/*")

foreach(_CHILD ${_SRC_CHILDREN})
    if(IS_DIRECTORY ${_CHILD})
        file(GLOB _CHILD_HEADERS "${_CHILD}/*.hpp")
	    file(GLOB _CHILD_SOURCES "${_CHILD}/*.cpp")
        list(APPEND MAIN_SOURCES ${_CHILD_HEADERS})
        list(APPEND MAIN_SOURCES ${_CHILD_SOURCES})
	endif()
endforeach()
###################################################

###################################################
# Create target and set properties
add_library(libmke2
    STATIC
        ${MAIN_HEADERS}
        ${MAIN_SOURCES}
)

target_include_directories(libmke2
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        "src"
)

add_executable(mke2 ${MAIN_EXECUTABLE})

target_link_libraries(mke2
    PUBLIC
        libmke2
)
###################################################
