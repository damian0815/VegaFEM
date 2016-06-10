cmake_minimum_required(VERSION 3.5)
project(sparseMatrixCmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-format-security")

set(SOURCE_FILES
    sparseMatrix.cpp
    sparseMatrix.h
    sparseMatrixIndexRemapper.cpp
    sparseMatrixIndexRemapper.h
    sparseMatrixMT.cpp
    sparseMatrixMT.h
    sparseMatrixOutline.cpp
    sparseMatrixOutline.h
    sparseSubMatrixLinkage.cpp
    sparseSubMatrixLinkage.h)

add_library(sparseMatrixCmake ${SOURCE_FILES})
