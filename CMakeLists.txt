cmake_minimum_required(VERSION 3.16)
project(parametron_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(parametron_ode src/ode.cpp)
target_compile_options(parametron_ode PRIVATE -O2 -Wall -Wextra -Wpedantic)

add_executable(parametron_logic src/logic.cpp)
target_compile_options(parametron_logic PRIVATE -O2 -Wall -Wextra -Wpedantic)

add_executable(parametron_tick src/tick.cpp)
target_compile_options(parametron_tick PRIVATE -O2 -Wall -Wextra -Wpedantic)

add_executable(phazeride_adc src/adc.cpp)
target_compile_options(phazeride_adc PRIVATE -O2 -Wall -Wextra -Wno-pedantic)

add_library(net_tick src/net_tick.cpp)
target_compile_options(net_tick PRIVATE -O2 -Wall -Wextra -Wpedantic)
target_include_directories(net_tick PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(parametron_tick PRIVATE net_tick)
target_link_libraries(phazeride_adc PRIVATE net_tick)
