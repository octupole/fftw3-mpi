# Policy
cmake_minimum_required( VERSION 3.0) 


# Project name
project(fftw3pp CXX C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake)
include(macros)
include(LookUp-GreatCMakeCookOff)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
option(PARALLEL "Enable Parallel compilation " ON)

if (PARALLEL)
   find_package(MPI)
   if (NOT MPI_CXX_FOUND)
 	  message(STATUS "========== OpenTRAJ Parallel Information ==========")
	  message(STATUS )
      message(STATUS "No MPI compiler found. Perhaps you want to set one explicitly?")
      message(STATUS "To override the default compiler, set the environment variable")
      message(STATUS "     export CXX=/full/path/to/mpic++")
      message(STATUS "in bash or ksh or sh.  In csh or tcsh, use:")
      message(STATUS "     setenv CXX /full/path/to/mpic++")
      message(STATUS "before running the cmake command.")
      message(STATUS )
   	  message(STATUS "==================================================")
   endif()
   IF ( MPI_CXX_FOUND )
   	  if(${CMAKE_MINOR_VERSION} GREATER_EQUAL 10)
   	  	if(NOT MPI_MPICXX_FOUND)
   	  		message(FATAL_ERROR "Your MPI implementation does not accept MPI c++ bindings. Connot compile.")
   	  	endif()
   	  elseif()
   	  	message(WARNING "Cannot figure out if your MPI implementation accept MPI c++ bindings. If it does not, it will not compile. Wait and see.")
   	  endif()
      INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
      if(DEFINED $ENV{CXX})
        message(STATUS "Found MPI compiler, but using the CXX compiler set in the CXX environment variable, $ENV{CXX}")
      else()
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        set(ENV{CXX} ${MPI_CXX_COMPILER})
      endif()
	  add_definitions(-DHAVE_MPI=1)
   ENDIF()
   include(FindOpenMP)
endif()
if(MPI_CXX_FOUND AND OpenMP_CXX_FLAGS)
	find_package(FFTW3 COMPONENTS DOUBLE MPI OPENMP)
	set(FFTW3_MPI_FOUND YES)
elseif(MPI_CXX_FOUND)
	find_package(FFTW3 COMPONENTS DOUBLE MPI)
	set(FFTW3_MPI_FOUND YES)
elseif(OpenMP_CXX_FLAGS)
	find_package(FFTW3 COMPONENTS DOUBLE OPENMP)
else()
	find_package(FFTW3 COMPONENTS DOUBLE)
endif()

if (OpenMP_CXX_FLAGS)
	add_definitions(-DHAVE_OPENMP=1)
endif()
if(FFTW3_MPI_FOUND)
	add_definitions(-DHAVE_FFTW_MPI=1)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
        message(FATAL_ERROR "Clang version must be at least 3.2!")
    endif()
    set(CMAKE_CXX_FLAGS "-std=c++11 -O3 ${OpenMP_CXX_FLAGS}")        ## Optimize	
    message(WARNING "Not tested, use at your won risk.")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
   set(CMAKE_CXX_FLAGS "-std=c++11 -ffast-math -fno-math-errno ${OpenMP_CXX_FLAGS}")  
   message(FATAL_ERROR "Cannot compile, AppleClang up to LLVM version 9.1.0 gives an internal error!")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0)
       message(FATAL_ERROR "Intel icc version must be >= 17. It is ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
   set(CMAKE_CXX_FLAGS "-std=c++11 ${OpenMP_CXX_FLAGS}")  

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.4)
        message(FATAL_ERROR "GCC version must be at least 5.4")
    endif()
   set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -fno-math-errno ${OpenMP_CXX_FLAGS}")  
endif()




if (FFTW3_DOUBLE_MPI_FOUND)
  include_directories (${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${FFTW3_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIRS})
  LINK_LIBRARIES(${FFTW3_LIBRARIES})
ELSE(FFTW3_FOUND)
  MESSAGE(WARNING "No fftw3 found - fftw3pp-mpi cannot be compiled")
endif (FFTW3_FOUND)

LONG_TODAY(BUILD_DATE)
  
add_subdirectory(fftw3pp)
add_subdirectory(comm)

if(FFTW3_FOUND)
	add_executable(test-fftw3 test-fftw3.cpp)
endif()


target_link_libraries(test-fftw3 fftw3pp comm)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "default install path" FORCE )
endif()
message( STATUS)
message( STATUS "========== OpenTRJ Build Information ===========")
message( STATUS "Current date/time .......... = ${BUILD_DATE}")
message( STATUS "CMAKE_SYSTEM ............... = ${CMAKE_SYSTEM}")
message( STATUS "================================================")
message( STATUS "CMAKE_BUILD_TYPE ........... = ${CMAKE_BUILD_TYPE}")
message( STATUS "CMAKE_INSTALL_PREFIX ....... = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "CMAKE_C_COMPILER ........... = ${CMAKE_C_COMPILER}")
message( STATUS "CMAKE_CXX_COMPILER ......... = ${CMAKE_CXX_COMPILER}")
if(MPI_CXX_FOUND)
message( STATUS "MPI_CXX_COMPILER ........... = ${MPI_CXX_COMPILER}")
message( STATUS "MPI_CXX_INCLUDE_PATH ....... = ${MPI_CXX_INCLUDE_PATH}")
message( STATUS "MPI_CXX_LIBRARIES .......... = ${MPI_CXX_LIBRARIES}")
message( STATUS "MPI_CXX_VERSION ............ = ${MPI_CXX_VERSION}")
message( STATUS "MPI_COMPILE_FLAGS .......... = ${CMAKE_CXX_FLAGS}")
endif()
if (OpenMP_CXX_FLAGS)
message(STATUS  "Compiled with OPENMP ....... = YES")
else()
message(STATUS  "Compiled with OPENMP ....... = NO ")
endif()
if(FFTW3_FOUND)
message( STATUS "FFTW3_INCLUDE_DIR .......... = ${FFTW3_INCLUDE_DIRS}")
message( STATUS "FFTW3_LIBRARIES ............ = ${FFTW3_LIBRARIES}")
if(FFTW3_DOUBLE_MPI_FOUND) 
message( STATUS "FFTW3_DOUBLE_MPI ............ = ${FFTW3_DOUBLE_MPI_LIBRARY}")
endif()
if(FFTW3_DOUBLE_OPENMP_FOUND) 
message( STATUS "FFTW3_DOUBLE_OPENMP ......... = ${FFTW3_DOUBLE_OPENMP_LIBRARY}")
endif()
if(FFTW3_DOUBLE_THREADS_FOUND) 
message( STATUS "FFTW3_DOUBLE_THREADS ........ = ${FFTW3_DOUBLE_THREADS_LIBRARY}")
endif()
else()
message( STATUS  "FFTW3            .......... = NOT FOUND")
endif()
message( STATUS)
message( STATUS "To override these options, add -D{OPTION_NAME}=... to the cmake command" )
message( STATUS "Particularly useful defines are for:")
message( STATUS )
message( STATUS "  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}	(where OpenTRJ executables will be installed)")
message( STATUS "  -DHINT_DIRS=/path/to/fftw3	(select path to installed FFTW3 routines)")
message( STATUS "  -DPARALLEL=OFF	(Do not compile with the default PARALLEL=ON)")
message( STATUS )

