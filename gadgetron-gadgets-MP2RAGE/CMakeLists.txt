cmake_minimum_required(VERSION 3.10...3.12 FATAL_ERROR)

project(GADGETRON-RMSB LANGUAGES CXX C)

find_package(Gadgetron REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(gadgets_MP2RAGE)