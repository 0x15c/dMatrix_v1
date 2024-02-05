# CMakeLists.txt in the root directory

cmake_minimum_required(VERSION 3.10)
project(dMatrix-v1)

# Add the include directory
include_directories(inc)

# Add subdirectories for source code
add_subdirectory(src)

# Add the main executable
add_executable(main main.c)

# Link the main executable with the my_math_library
target_link_libraries(main PRIVATE dMatrix)
