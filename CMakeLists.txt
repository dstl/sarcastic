cmake_minimum_required(VERSION 3.11.4)
project(sarcastic)
set(CMAKE_FIND_FRAMEWORK LAST)
set(NO_CMAKE_SYSTEM_PATH)
set(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS "-std=c++11")

# Make sure git exists and use it to get version information
#
find_package(Git)
execute_process(COMMAND bash "-c" "${GIT_EXECUTABLE} describe --dirty | sed -e 's/^v//' -e 's/^V//' -e 's/g//'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE res
		OUTPUT_VARIABLE FULL_VERSION
		ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} describe --abbrev=0 | sed -e 's/^v//' -e 's/^V//'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	SHORT_VERSION
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} rev-list master | wc -l | awk '{print $1}'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	REVISION
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} show --format=\"%ci\" | head -1"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	VERSION_DATE
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)

# Configure version header files
#
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/Version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/Version.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SartraceVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SartraceVersion.h)

# Check for required packages
#
find_package(gdal 2.3.0 REQUIRED)
if(${GDAL_FOUND})
    include_directories(SARTrace, ${GDAL_INCLUDE_DIR})
    add_definitions(-Dsarclib_HAVE_GDAL) 
endif()



#  Set file specific add_definitions
#
set_property( SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/fft_c8.c PROPERTY COMPILE_DEFINITIONS TRUNCATE )


#The file(GLOB...) allows for wildcard additions:
file(GLOB SARCLIBSOURCES "src/sarclib/*.c")
file(GLOB SARTraceSOURCES "src/SARTrace/*.c*")
 
#Generate the library from the sources
add_library(sarclib STATIC ${SARCLIBSOURCES})
add_executable(SARTrace ${SARTraceSOURCES})
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic)

# find_package(opencl REQUIRED)
# if(${OpenCL_FOUND})
# 	include_directories( ${OPENCL_INCLUDE_DIR} )
# 	target_include_directories(SARTrace PRIVATE ${OpenCL_INCLUDE_DIRS})
# 	message(" OPENCL include dir: ${OpenCL_INCLUDE_DIRS}")
# endif()

#Set the location for library installation
install(TARGETS sarclib DESTINATION src/lib)

#target_link_libraries(sarclib gdal)