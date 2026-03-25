set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++)

cmake_minimum_required(VERSION 3.18)

project(nicole)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(odepack_cpp)
add_subdirectory(nicole)
