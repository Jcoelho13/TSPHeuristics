cmake_minimum_required(VERSION 3.24)
project(DATP2)

set(CMAKE_CXX_STANDARD 17)

add_executable(DATP2 main.cpp DataStructures/graph.cpp DataStructures/graph.h DataStructures/vertexEdge.cpp DataStructures/vertexEdge.h DataStructures/MutablePriorityQueue.h data/data.h DataStructures/Node.cpp DataStructures/Node.h)
