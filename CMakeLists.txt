cmake_minimum_required(VERSION 3.11.4)
project(sarcastic)
set(NO_CMAKE_SYSTEM_PATH)
set(CMAKE_BUILD_TYPE Release)
find_package(gdal 2.3.0 REQUIRED)
message("version: ${GDAL_VERSION}")

if(${GDAL_FOUND})
    message("****GDAL inc dir is at ${GDAL_INCLUDE_DIR}")
    include_directories(include, ${GDAL_INCLUDE_DIR})
    add_definitions(-Dsarclib_HAVE_GDAL) 
endif()
 
#The file(GLOB...) allows for wildcard additions:
file(GLOB SARCLIBSOURCES "src/sarclib/*.c")
 
#Generate the library from the sources
add_library(sarclib STATIC ${SARCLIBSOURCES})
 
#Set the location for library installation
install(TARGETS sarclib DESTINATION src/lib)

#target_link_libraries(sarclib gdal)