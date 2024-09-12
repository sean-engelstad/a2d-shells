# Specify cmake version needed
cmake_minimum_required(VERSION 3.12)

project("A2D_SHELLS"
  VERSION 0.0.1
  DESCRIPTION "MiniApp for A2D Implementation of Geometric nonlinear shells"
  HOMEPAGE_URL "https://github.com/sean-engelstad/a2d-shells"
  LANGUAGES CXX)

# Export the compile commands that clangd LS uses
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# C++ 11 is required for TACS / A2D-SHELLS
set(CMAKE_CXX_STANDARD 11)

set(A2D_DIR $ENV{HOME}/installs/a2d)

# If in debug mode, set the preprocessor definition
if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DASHELLS_DEBUG_MODE)
endif()

# Dependencies
include(cmake/utils.cmake)
# check_path_and_inform(FATAL_ERROR A2D ${A2D_DIR} A2D_DIR)
find_package(A2D REQUIRED PATHS ${A2D_DIR})

# Set warning flags
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|GNU")
  add_compile_options(
    -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
  )
endif()

# Copy over inputs
# file(COPY ${PROJECT_SOURCE_DIR}/input DESTINATION ${PROJECT_BINARY_DIR})

# === Install A2D_SHELLS ===

# Change the default value of CMAKE_INSTALL_PREFIX
set(A2D_SHELLS_INSTALL_DIR $ENV{HOME}/installs/a2d-shells)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${A2D_SHELLS_INSTALL_DIR} CACHE PATH "The default\
  path to install xcgd to" FORCE)
  message("The directory to install a2d-shells is set to default:\
  ${A2D_SHELLS_INSTALL_DIR}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Interface makes the target header-only that does not need to be compiled
add_library(${PROJECT_NAME} INTERFACE)

# Pass along c++11 requirements to projects using this library
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Dependencies
# target_include_directories(${PROJECT_NAME} INTERFACE
#  ${XCGD_ALGOIM_DIR}/algoim)

# if (A2D_SHELLS_INSTALL_LIBRARY)
install(
TARGETS ${PROJECT_NAME}
EXPORT ${PROJECT_NAME}_Targets
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
VERSION ${PROJECT_VERSION}
COMPATIBILITY SameMajorVersion)

message(${PROJECT_NAME}Config.cmake.in)

configure_package_config_file(
"${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
INSTALL_DESTINATION
${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# install cmake files
install(
EXPORT ${PROJECT_NAME}_Targets
FILE "${PROJECT_NAME}Targets.cmake"
NAMESPACE ${PROJECT_NAME}::
DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)
install(
FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
"${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# Install headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)