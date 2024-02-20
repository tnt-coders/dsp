cmake_minimum_required(VERSION 3.26.0)

# Define the project
project(dsp LANGUAGES CXX)

# Include required modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(tnt/project)

# Set project attributes
tnt_project__set_version_from_git()
tnt_project__set_namespace(tnt)

# Enable testing for the project
# Note: must be in top level CMakeLists.txt
enable_testing()

add_subdirectory(docs)
add_subdirectory(src)
add_subdirectory(test)

# Install project targets
tnt_project__install(TARGETS ${PROJECT_NAME})
