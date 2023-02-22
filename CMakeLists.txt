cmake_minimum_required(VERSION 3.2)

project(WeldFormSPH)
enable_language(Fortran)
                            
                            
INCLUDE      (FindOpenMP)

if(OPENMP_FOUND)
    ADD_DEFINITIONS (-DUSE_OMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
        SET (MISSING "${MISSING} OpenMP")
endif(OPENMP_FOUND)

#if GNU, 
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall"})
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()

file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(WeldFormSPH ${sources})
add_executable (openmp examples/omptest.f90)
add_executable (mechtest mechtest.f90
                          src/Domain.f90
                          src/Functions.f90
                          src/Integrate.f90
                          src/Kernel.f90
                          src/Mechanical.f90
                          src/Neighbor.f90
                          src/ParticleData.f90
                          src/Precision.f90
                          src/SolverFraser.f90
                          )




# project(WeldFormGPUSPH Fortran CUDA)
# find_package(CUDAToolkit REQUIRED)
# target_link_libraries(program CUDA::cusolver CUDA::cudart)

# set_target_properties(program PROPERTIES CUDA_ARCHITECTURES 70) # <---This doesn't seem to work either
# set_property(TARGET program PROPERTY CUDA_ARCHITECTURES 70) # <---This doesn't seem to work either

# set_target_properties(program PROPERTIES LINK_FLAGS "-Mcuda") # -gpu=cc70") # <---Only adding these flags, it compiles for the desired architectures. Also, adding the "-Mcuda" is needed.
# string(CONCAT FORTRAN_FLAGS "-cpp "
                            # "-lineinfo "
                            # "-O3 "
                            # #"-gpu=cc70 " # <-- This is also needed
                            # "-Minfo=accel "
                            # "-Mcuda " # <-- This is also needed
                            # "-g "
                            # )