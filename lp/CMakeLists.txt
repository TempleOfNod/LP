cmake_minimum_required(VERSION 3.1)
project(lp)

add_executable(${PROJECT_NAME} lp.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)