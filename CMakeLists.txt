cmake_minimum_required(VERSION 3.14)
project(ccmlina)

# ---------------------#
# --- Compile Flags ---#
# ---------------------#
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
enable_language(C CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -Wno-sign-compare -Wno-unused-parameter -fopenmp")

    # Build-type specific flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

    # Static libs
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

message(STATUS "Building with the following extra flags: ${CMAKE_CXX_FLAGS}")

# -------------#
# --- Files ---#
# -------------#
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SOURCES "src/*.cpp*")

# ---------------#
# --- Library ---#
# ---------------#
add_library(${PROJECT_NAME} SHARED ${SOURCES})

# ---------------#
# --- Testing ---#
# ---------------#
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
file(GLOB TESTSOURCES "tests/*.cpp*")
add_executable(
    mtest
    ${TESTSOURCES}
)
target_link_libraries(
    mtest
    GTest::gtest_main
    ${PROJECT_NAME}
)

include(GoogleTest)
gtest_discover_tests(mtest)