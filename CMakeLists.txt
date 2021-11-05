cmake_minimum_required(VERSION 3.10)
project(example)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# polyscope
add_subdirectory("external/polyscope")

# fabsim
add_subdirectory("external/fabsim")

# optim
option(OPTIM_USE_CHOLMOD  "Use Cholmod linear solver instead of Eigen's built-in one"  OFF)
option(OPTIM_USE_PARDISO  "Compile Optim with Intel Pardiso"  OFF)
add_subdirectory("external/optim")

# c++ flags
set(CMAKE_CXX_STANDARD 17)

# Project files
add_executable(membrane_example src/membrane.cpp)
add_executable(rod_example src/rod.cpp)

target_link_libraries(membrane_example polyscope fabsim optim)
target_link_libraries(rod_example polyscope fabsim optim)
