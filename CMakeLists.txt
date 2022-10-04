cmake_minimum_required(VERSION 3.2)

project(inverse_kinematics)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(inverse_kin)
add_subdirectory(matrix)

add_executable(inv_kin_exec main.cpp)
target_link_libraries(inv_kin_exec matrix inverse_kin)
