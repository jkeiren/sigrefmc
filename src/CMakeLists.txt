cmake_minimum_required(VERSION 2.6)
project(sigref C CXX)

# Check for GMP
include(FindGMP)
if(GMP_FOUND)
  message(STATUS "GMP Libraries: " ${GMP_LIB} " " ${GMPXX_LIB} " with include path " ${GMP_INCLUDE_DIR} " " ${GMPXX_INCLUDE_DIR})
  include_directories(${GMP_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Cannot find GMP library")
endif()

# Check for TinyXML
include(FindTinyXML)
if(TINYXML_FOUND)
  message(STATUS "TinyXML Libraries: " ${TINYXML_LIB} " with include path " ${TINYXML_INCLUDE_DIR})
  include_directories(${TINYXML_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Cannot find TinyXML library")
endif()

include_directories(.)

set(SOURCES
    bisimulation.hpp
    bisim_ctmc.cpp
    bisim_lts.cpp
    bisim_lts2.cpp
    bisim_imc.cpp
    getrss.h
    getrss.c
    blocks.h
    blocks.c
    inert.h
    inert.c
    parse_bdd.hpp
    parse_bdd.cpp
    parse_xml.hpp
    parse_xml.cpp
    refine.h
    systems.hpp
    sigref.h
    sigref.cpp
    sigref_util.hpp
    sigref_util.cpp
    quotient.hpp
    quotient.cpp
    writer.hpp
    writer.cpp
    )

add_executable(sigrefmc ${SOURCES} refine_sl.c)
target_link_libraries(sigrefmc sylvan ${TINYXML_LIB} ${GMP_LIB})

add_executable(sigrefmc_ht ${SOURCES} refine_ht.c)
target_link_libraries(sigrefmc_ht sylvan ${TINYXML_LIB} ${GMP_LIB})

include(CheckIncludeFiles)
check_include_files("gperftools/profiler.h" HAVE_PROFILER)

find_package(Boost REQUIRED)


if(HAVE_PROFILER)
    set_target_properties(sigrefmc PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(sigrefmc profiler)

    set_target_properties(sigrefmc_ht PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(sigrefmc_ht profiler)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # add argp library for OSX
    target_link_libraries(sigrefmc argp)
    target_link_libraries(sigrefmc_ht argp)
endif()
