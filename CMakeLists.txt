cmake_minimum_required(VERSION 3.0)
project(mke2 VERSION 0.1.0 LANGUAGES CXX)

set(ENABLE_TESTS "ON")
set(ENABLE_OPENMP "OFF")

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(LIBPROXY "mkeproxyd")
    set(EXECUTABLE "mked")
else()
    set(LIBPROXY "mkeproxy")
    set(EXECUTABLE "mke")
endif()

###################################################
# Include additional cmake packages
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
###################################################

set(CMAKE_CXX_STANDARD 17)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

###################################################
# Configure conan
if (ENABLE_TESTS)
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/conan.cmake)
        message(STATUS
            "Downloading conan.cmake from https://github.com/conan-io/cmake-conan"
        )
        file(DOWNLOAD
            "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake"
        )
    endif()
    include(${CMAKE_BINARY_DIR}/conan.cmake)
    conan_cmake_run(
        REQUIRES
            gtest/1.8.1@bincrafters/stable
        BASIC_SETUP CMAKE_TARGETS
        BUILD missing
        OPTIONS
    )
endif()
###################################################

###################################################
# Declare dependencies
if (ENABLE_OPENMP)
    find_package(OpenMP)
endif()

if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /analyze")
endif()

# add_subdirectory(gpu_solver)
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
add_library(${LIBPROXY}
    STATIC
        ${MAIN_HEADERS}
        ${MAIN_SOURCES}
)

target_include_directories(${LIBPROXY}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        "src"
)

target_link_libraries(${LIBPROXY}
    PRIVATE
        # gpu_lib
)

add_executable(${EXECUTABLE} ${MAIN_EXECUTABLE})

target_link_libraries(${EXECUTABLE}
    PUBLIC
        ${LIBPROXY}
)
###################################################

if(ENABLE_TESTS)
    add_subdirectory(test)
endif()