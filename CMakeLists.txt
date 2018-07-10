cmake_minimum_required(VERSION 3.11.4)
project(sarcastic)
set(CMAKE_FIND_FRAMEWORK LAST)
set(NO_CMAKE_SYSTEM_PATH)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/SarcasticVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/SarcasticVersion.h)

# Check for required packages
#
find_package(gdal 2.3.0 REQUIRED)
if(${GDAL_FOUND})
    include_directories(SARTrace, ${GDAL_INCLUDE_DIR})
    add_definitions(-Dsarclib_HAVE_GDAL) 
endif()

find_package(opencl REQUIRED)
if(${OpenCL_FOUND})
	include_directories( ${OPENCL_INCLUDE_DIR} )
	message(" OPENCL include dir: ${OpenCL_INCLUDE_DIRS}")
	message(" OpenCL_LIBRARIES: ${OpenCL_LIBRARIES}")
endif()

find_package(expat REQUIRED)
if(${EXPAT_FOUND})
	message(" EXPAT_LIBRARIES: ${EXPAT_LIBRARIES}")
endif()

find_package(fftw REQUIRED)
if(${FFTW_FOUND})
	include_directories( ${FFTW_INCLUDE_DIR} )
	set(CMAKE_EXE_LINKER_FLAGS "-lfftw3f" CACHE INTERNAL "" FORCE)
endif()

#  Set file specific add_definitions
#
set_property( SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/fft_c8.c PROPERTY COMPILE_DEFINITIONS TRUNCATE )


#The file(GLOB...) allows for wildcard additions:
file(GLOB SARCLIBSOURCES "src/sarclib/*.c")
 
add_library(sarclib STATIC ${SARCLIBSOURCES})


#Generate the library from the sources
add_executable(SARTrace
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/AABB.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/clipToAABB.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/rayTrace.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/kdTreeNode.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/splitCandidate.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinystr.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxml.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlerror.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlparser.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/buildTree.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/colladainterface.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/TriangleMesh.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRays.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/accelerateTriangles.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRopesAndBoxes.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/buildKernel.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/reflect.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/getSARTraceUserInput.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SARTrace.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SARTracebanner.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/traceThreadCore.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/readMaterialFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/matrixMultiplication.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/OpenCLUtils.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/ecef2SceneCoords.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/pmath_rng.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/boxMullerRandom.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/ranf.c
	)

target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/bircs)
target_include_directories(SARTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile)
target_link_libraries(SARTrace sarclib)
target_link_libraries(sarclib ${FFTW_LIBRARIES})
target_link_libraries(sarclib ${EXPAT_LIBRARIES})
target_link_libraries(SARTrace ${OpenCL_LIBRARIES})



#Set the location for library installation
install(TARGETS sarclib DESTINATION src/lib)

#target_link_libraries(sarclib gdal)