# ChemLarp — cross-platform build (Windows: MSVC or clang-cl; macOS: Apple clang; Linux: GCC/Clang).
#
#   cmake -S . -B build
#   cmake --build build --config Release --parallel
#   ctest --test-dir build -C Release --output-on-failure
cmake_minimum_required(VERSION 3.20)
project(chemlarp_core VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Put the DLL / dylib, CLI and tests side by side (Windows finds DLLs next to the .exe).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(CHEMLARP_BUILD_TESTS "Build the test suite" ON)
option(CHEMLARP_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

add_library(chemlarp_warnings INTERFACE)
if(MSVC)
  # /utf-8: sources and CLI output contain UTF-8; /permissive-: standards conformance.
  target_compile_options(chemlarp_warnings INTERFACE /W4 /utf-8 /permissive- /Zc:__cplusplus /EHsc
                         /wd4244 /wd4267 /wd4305)  # implicit numeric narrowing is intentional in hot loops
  target_compile_definitions(chemlarp_warnings INTERFACE _CRT_SECURE_NO_WARNINGS NOMINMAX)
  if(CHEMLARP_WARNINGS_AS_ERRORS)
    target_compile_options(chemlarp_warnings INTERFACE /WX)
  endif()
else()
  target_compile_options(chemlarp_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow)
  if(CHEMLARP_WARNINGS_AS_ERRORS)
    target_compile_options(chemlarp_warnings INTERFACE -Werror)
  endif()
endif()

file(GLOB CHEMLARP_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Static library: linked into the CLI and tests.
find_package(Threads REQUIRED)
add_library(chemlarp_static STATIC ${CHEMLARP_SOURCES})
target_include_directories(chemlarp_static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_compile_definitions(chemlarp_static PUBLIC CHEMLARP_STATIC)
target_link_libraries(chemlarp_static PRIVATE chemlarp_warnings PUBLIC Threads::Threads)

# Shared library: the C ABI for Python (ctypes), Node and other FFI consumers.
#   Windows: chemlarp.dll (+ chemlarp.lib import library)   macOS: libchemlarp.dylib   Linux: libchemlarp.so
add_library(chemlarp SHARED ${CHEMLARP_SOURCES})
target_include_directories(chemlarp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_compile_definitions(chemlarp PRIVATE CHEMLARP_BUILD_DLL)
target_link_libraries(chemlarp PRIVATE chemlarp_warnings Threads::Threads)
set_target_properties(chemlarp PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
if(MINGW)
  # Self-contained DLL: no libstdc++/libgcc DLLs needed next to it.
  target_link_options(chemlarp PRIVATE -static-libgcc -static-libstdc++ -static)
endif()

add_executable(chemlarp_cli src/cli/main.cpp)
set_target_properties(chemlarp_cli PROPERTIES OUTPUT_NAME chemlarp)
target_link_libraries(chemlarp_cli PRIVATE chemlarp_static chemlarp_warnings)

if(CHEMLARP_BUILD_TESTS)
  enable_testing()
  file(GLOB CHEMLARP_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
  add_executable(chemlarp_tests ${CHEMLARP_TESTS})
  target_link_libraries(chemlarp_tests PRIVATE chemlarp_static chemlarp_warnings)
  target_compile_definitions(chemlarp_tests PRIVATE CHEMLARP_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
  add_test(NAME chemlarp_tests COMMAND chemlarp_tests WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

include(GNUInstallDirs)
install(TARGETS chemlarp chemlarp_cli
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/chemlarp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/chemlarp)
