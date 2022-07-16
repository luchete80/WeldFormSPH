cmake_minimum_required(VERSION 3.2)

project(WeldFormSPH)
enable_language(Fortran)

INCLUDE      (FindOpenMP)

file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(WeldFormSPH ${sources})
