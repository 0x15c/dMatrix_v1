# CMakeLists.txt in the src directory

# Add source files
add_library(dMatrix STATIC dComplex.c matrix.c)

# Include the headers from the inc directory
target_include_directories(dMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
