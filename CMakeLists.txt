cmake_minimum_required(VERSION 3.10)
project(example)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_subdirectory("external/fabsim")
add_subdirectory("external/optim")
add_subdirectory("external/polyscope")

# c++ flags
set(CMAKE_CXX_STANDARD 17)

# Project files
add_executable(composite_example src/composite.cpp)
add_executable(membrane_example src/membrane.cpp)
add_executable(rod_example src/rod.cpp)

target_link_libraries(composite_example fabsim optim polyscope)
target_link_libraries(membrane_example fabsim optim polyscope)
target_link_libraries(rod_example fabsim optim polyscope)
