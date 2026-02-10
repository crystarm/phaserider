cmake_minimum_required(VERSION 3.16)
project(parametron_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(parametron_ode src/ode.cpp)
target_compile_options(parametron_ode PRIVATE -O2 -Wall -Wextra -Wpedantic)
