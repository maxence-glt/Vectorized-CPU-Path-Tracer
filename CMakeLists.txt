cmake_minimum_required(VERSION 3.20)
project(raytracer LANGUAGES CXX)

option(ENABLE_SANITIZERS "Enable ASan/UBSan/LSan in Debug builds" ON)
option(BUILD_TESTS "Enable unit test (Google Test)" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CORE_SRCS
    src/worlds/manyballs.cpp
    src/util/error.cpp
    src/util/log.cpp
    src/util/profiler.cpp
    src/util/transform.cpp
    src/render/render.cpp
    src/raytracer.cpp
    src/util/color.cpp
)

add_library(raytracer_core STATIC ${CORE_SRCS})

# ---- OpenImageIO ----
set(OpenImageIO_ROOT "${CMAKE_SOURCE_DIR}/import/OpenImageIO-3.1.9.0/install")
list(PREPEND CMAKE_PREFIX_PATH "${OpenImageIO_ROOT}")

find_package(OpenImageIO CONFIG REQUIRED)
target_link_libraries(raytracer_core PUBLIC OpenImageIO::OpenImageIO)

target_include_directories(raytracer_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# release flags
target_compile_options(raytracer_core PRIVATE
  -fmacro-prefix-map=${CMAKE_SOURCE_DIR}/=
  -fmacro-prefix-map=${CMAKE_BINARY_DIR}/=
)

target_compile_options(raytracer_core PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-march=native>
)

# debug flags
target_compile_definitions(raytracer_core PRIVATE
    $<$<CONFIG:Debug>:DEBUG_BUILD>
)

target_compile_options(raytracer_core PRIVATE
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
    $<$<CONFIG:Debug>:-fno-optimize-sibling-calls>
    $<$<CONFIG:Debug>:-Wno-character-conversion>
)

# main executable
add_executable(raytracer src/main.cpp)
target_link_libraries(raytracer PRIVATE raytracer_core)

set_target_properties(raytracer PROPERTIES
    OUTPUT_NAME "raytracer"
    OUTPUT_NAME_DEBUG "raytracer_debug"
)

if(ENABLE_SANITIZERS)
    add_library(sanitizers INTERFACE)
    target_compile_options(sanitizers INTERFACE
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,leak>
    $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
  )
    target_link_options(sanitizers INTERFACE
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,leak>
  )

    # apply to EVERYTHING that links into a final executable
    target_link_libraries(raytracer_core PUBLIC sanitizers)
    target_link_libraries(raytracer PRIVATE sanitizers)

    if(TARGET tests)
        target_link_libraries(tests PRIVATE sanitizers)
    endif()
endif()

# critical, stable dependencies
set(OpenImageIO_ROOT "${CMAKE_SOURCE_DIR}/import/OpenImageIO-3.1.9.0/install")
list(PREPEND CMAKE_PREFIX_PATH "${OpenImageIO_ROOT}")
find_package(OpenImageIO CONFIG REQUIRED)

target_link_libraries(raytracer_core PUBLIC OpenImageIO::OpenImageIO)

# development/testing dependencies
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest CONFIG REQUIRED)

    file(GLOB_RECURSE TEST_SRCS CONFIGURE_DEPENDS test/*.cpp)

    add_executable(tests EXCLUDE_FROM_ALL ${TEST_SRCS})
    target_link_libraries(tests PRIVATE raytracer_core GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(tests)
endif()

target_compile_options(raytracer PRIVATE
  -fmacro-prefix-map=${CMAKE_SOURCE_DIR}/=
  -fmacro-prefix-map=${CMAKE_BINARY_DIR}/=
)

