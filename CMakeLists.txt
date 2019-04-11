cmake_minimum_required(VERSION 3.11.4)
project(sarcastic)
# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ELSE()
	find_package (Threads)
ENDIF()

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_FIND_FRAMEWORK LAST)
set(NO_CMAKE_SYSTEM_PATH)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -D_GNU_SOURCE")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

add_definitions(-DKERNELDIR="${CMAKE_SOURCE_DIR}/src/kernels") 
add_definitions(-DMATERIALPROPS="${CMAKE_SOURCE_DIR}/scripts/MaterialProperties.txt") 

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
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} log -1 --pretty=format:%h"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	COMMIT
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} describe | sed -e 's/^v//' -e 's/^V//' | awk -F. '{print $1}'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	MAJOR
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} describe | sed -e 's/^v//' -e 's/^V//' | awk -F. '{print $2}' | awk -F- '{print $1}'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	MINOR
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND	bash "-c" "${GIT_EXECUTABLE} describe | sed -e 's/^v//' -e 's/^V//' | awk -F. '{print $2}' | awk -F- '{print $2}'"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		RESULT_VARIABLE	res
		OUTPUT_VARIABLE	PATCH
		ERROR_QUIET	OUTPUT_STRIP_TRAILING_WHITESPACE)

message("FULL_VERSION is ${FULL_VERSION}")
message("SHORT_VERSION is ${SHORT_VERSION}")
message("REVISION is ${REVISION}")
message("MAJOR is ${MAJOR}")
message("MINOR is ${MINOR}")
message("PATCH is ${PATCH}")

# Configure version header files
#
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/Version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/Version.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SartraceVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/SARTrace/SartraceVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/SarcasticVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/SarcasticVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/BircsVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/BircsVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/CTTVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/CTTVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/materialise/materialise_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/materialise/materialise_version.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cphdInfo/readCPHD_Version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/cphdInfo/readCPHD_Version.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cphdShell/cphdShell_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/cphdShell/cphdShell_version.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/tdpoclVersion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/tdpoclVersion.h)

# Check for required packages
#
find_package(GDAL REQUIRED)
if(${GDAL_FOUND})
    include_directories(SARTrace, ${GDAL_INCLUDE_DIR})
    add_definitions(-Dsarclib_HAVE_GDAL) 
endif()

find_package(OpenCL REQUIRED)
if(${OpenCL_FOUND})
	include_directories( ${OPENCL_INCLUDE_DIR} )
endif()

find_package(EXPAT REQUIRED)
if(${EXPAT_FOUND})
	link_libraries(${EXPAT_LIBRARIES})
endif()

find_package(FFTW REQUIRED)
if(${FFTW_FOUND})
	include_directories( ${FFTW_INCLUDE_DIR} )
	set(CMAKE_EXE_LINKER_FLAGS "-L/usr/local/lib -lfftw3f -lfftw3f_threads" CACHE INTERNAL "" FORCE)
endif()

find_package(Boost REQUIRED COMPONENTS filesystem)
if(${BOOST_FOUND})
	include_directories(${BOOST_INCLUDE_DIR})
endif()

find_package(Readline REQUIRED)
if(${READLINE_FOUND})
    add_definitions(-DHAVE_READLINE=1)
    link_libraries(${READLINE_LIBRARIES})
endif()

find_package(CGAL REQUIRED)
if(${CGAL_FOUND})
	include_directories(${CGAL_INCLUDE_DIR})
endif()

find_package (Threads)

#  Set file specific add_definitions
#
set_property( SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarclib/fft_c8.c PROPERTY COMPILE_DEFINITIONS TRUNCATE )


#The file(GLOB...) allows for wildcard additions:
file(GLOB SARCLIBSOURCES "src/sarclib/*.c")
add_library(sarclib STATIC ${SARCLIBSOURCES})
include_directories(src)

#Generate SARTrace executable
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
target_include_directories(SARTrace PRIVATE ${OpenCL_INCLUDE_DIR})
target_link_libraries(SARTrace sarclib)
target_link_libraries(SARTrace ${FFTW_LIBRARIES})
target_link_libraries(SARTrace ${OpenCL_LIBRARIES})
target_link_libraries(SARTrace stdc++)
target_link_libraries(SARTrace m)

# Generate bircs executable
add_executable(bircs
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRays.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/rayTrace.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlparser.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinystr.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxml.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlerror.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/colladainterface.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/splitCandidate.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/AABB.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/clipToAABB.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/kdTreeNode.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/buildTree.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/TriangleMesh.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/tryReadFile.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/accelerateTriangles.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/buildKernel.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRopesAndBoxes.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/cpuPOField.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/main.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/threadCore.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/shadowRays.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/reflect.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/readMaterialFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/bircs/bircsBanner.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/TxPowerPerRay.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/matrixMultiplication.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/ranf.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/boxMullerRandom.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/pmath_rng.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/RCS.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/OpenCLUtils.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/printProgress.c
)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/bircs)
target_include_directories(bircs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile)
target_include_directories(bircs PRIVATE ${OpenCL_INCLUDE_DIR})
target_link_libraries(bircs sarclib)
target_link_libraries(bircs ${FFTW_LIBRARIES})
target_link_libraries(bircs ${OpenCL_LIBRARIES})
target_link_libraries(bircs ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bircs stdc++)
target_link_libraries(bircs m)

# Generate colladaToPlyFile executable
add_executable(colladaToPlyFile 
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/colladainterface.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/main.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/banner.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinystr.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxml.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlerror.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlparser.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/TriangleMesh.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/readMaterialFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/matrixMultiplication.c
)
target_include_directories(colladaToPlyFile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_include_directories(colladaToPlyFile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(colladaToPlyFile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile)
target_link_libraries(colladaToPlyFile sarclib)
target_link_libraries(colladaToPlyFile stdc++)
target_link_libraries(colladaToPlyFile m)

# Generate materialise executable
add_executable(materialise
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/colladainterface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/materialise/main.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinystr.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxml.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlerror.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlparser.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/TriangleMesh.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/readMaterialFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/boxMullerRandom.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/matrixMultiplication.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/pmath_rng.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/printProgress.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/ranf.c
)
target_include_directories(materialise PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(materialise PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile)
target_include_directories(materialise PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_link_libraries(materialise sarclib)
target_link_libraries(materialise ${Boost_LIBRARIES})
target_link_libraries(materialise boost_filesystem)
target_link_libraries(materialise ${GMP_LIBRARIES})
target_link_libraries(materialise ${CGAL_LIBRARIES})
target_link_libraries(materialise stdc++)
target_link_libraries(materialise m)

# Generate sarcastic executable
add_executable(sarcastic
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/AABB.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/banner.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRays.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/buildTree.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/clipToAABB.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile/colladainterface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/kdTreeNode.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/rayTrace.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/reflect.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree/splitCandidate.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinystr.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxml.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlerror.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/tinyxml/tinyxmlparser.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/TriangleMesh.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/getUserInput.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/shadowRays.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/accelerateTriangles.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/cpuPOField.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/buildRopesAndBoxes.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/sarcastic2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/tryReadFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/threadCore.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/readMaterialFile.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/boxMullerRandom.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/printProgress.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/sceneExtent.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/ecef2SceneCoords.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/matrixMultiplication.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/ranf.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/RCS.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/TxPowerPerRay.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/pmath_rng.c
)
target_include_directories(sarcastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic)
target_include_directories(sarcastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(sarcastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fastKdTree)
target_include_directories(sarcastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/colladaToPlyFile)
target_include_directories(sarcastic PRIVATE ${OpenCL_INCLUDE_DIR})
target_link_libraries(sarcastic sarclib)
target_link_libraries(sarcastic ${FFTW_LIBRARIES})
target_link_libraries(sarcastic ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(sarcastic stdc++)
target_link_libraries(sarcastic m)

# Generate cphdInfo executable
add_executable(cphdInfo
	${CMAKE_CURRENT_SOURCE_DIR}/src/cphdInfo/main.cpp
)
target_link_libraries(cphdInfo sarclib)
target_link_libraries(cphdInfo stdc++)
target_link_libraries(cphdInfo m)

# Generate cphdShell executable
add_executable(cphdShell
	${CMAKE_CURRENT_SOURCE_DIR}/src/cphdShell/main.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/cphdShell/ecef2SceneCoords.cpp
)
target_link_libraries(cphdShell sarclib)
target_link_libraries(cphdShell stdc++)
target_link_libraries(cphdShell m)

# Generate tdpocl executable
add_executable(tdpocl
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/tdp.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/getUserInput.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/createSurface.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/topotrimSetup.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/banner.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/tdpcore.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/getMemorySize.c 
	${CMAKE_CURRENT_SOURCE_DIR}/src/tdpocl/oclSRPRanges.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/Common/OpenCLUtils.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic/TxPowerPerRay.c
)
target_include_directories(tdpocl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Common)
target_include_directories(tdpocl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/sarcastic)
target_include_directories(tdpocl PRIVATE ${OpenCL_INCLUDE_DIR})
target_link_libraries(tdpocl sarclib)
target_link_libraries(tdpocl ${FFTW_LIBRARIES})
target_link_libraries(tdpocl ${OpenCL_LIBRARIES})
target_link_libraries(tdpocl m)

# Install the executables
install(TARGETS 
	SARTrace 
	bircs 
	colladaToPlyFile 
	cphdInfo 
	cphdShell 
	materialise 
	sarcastic 
	tdpocl
	DESTINATION bin
)

set(CPACK_PACKAGE_NAME "sarcastic")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PATCH}")
set(CPACK_SOURCE_IGNORE_FILES "/build;/.git/;/bin/;/.DS_Store;/.gitignore;/*.pyc;~$;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

#target_link_libraries(sarclib gdal)
