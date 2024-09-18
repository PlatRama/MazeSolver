cmake_minimum_required(VERSION 3.28)
project(MazeSolver)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")

add_executable(MazeSolver main.cpp
        Maze.cpp
        StatsUtil.cpp
)
