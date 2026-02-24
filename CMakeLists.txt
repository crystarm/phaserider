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

add_library(parametron_phys src/parametron_phys.cpp)
target_compile_options(parametron_phys PRIVATE -O2 -Wall -Wextra -Wpedantic)
target_include_directories(parametron_phys PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(parametron_phys_adder src/parametron_phys_adder.cpp)
target_compile_options(parametron_phys_adder PRIVATE -O2 -Wall -Wextra -Wpedantic)
target_include_directories(parametron_phys_adder PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(parametron_phys_adder PRIVATE parametron_phys)

add_library(circuit_phys src/circuit_phys.cpp)
target_compile_options(circuit_phys PRIVATE -O2 -Wall -Wextra -Wpedantic)
target_include_directories(circuit_phys PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(circuit_tick src/circuit_tick.cpp)
target_compile_options(circuit_tick PRIVATE -O2 -Wall -Wextra -Wpedantic)
target_include_directories(circuit_tick PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(parametron_tick PRIVATE circuit_tick)
target_link_libraries(phazeride_adc PRIVATE circuit_tick)
