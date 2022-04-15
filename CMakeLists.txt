cmake_minimum_required(VERSION 3.5)
project(MDMTraceExample)

set(CMAKE_CXX_STANDARD 17)

enable_language (Fortran)
set (CMAKE_Fortran_FLAGS "-fPIC -finit-local-zero -falign-commons -fno-automatic")
set (CMAKE_CXX_FLAGS "-fPIC")

# build libraray libJsonCpp
add_library(JsonCpp SHARED src/jsoncpp.cpp)
target_include_directories(JsonCpp PUBLIC ${PROJECT_SOURCE_DIR}/include/json)

# build library libMDMTraceFit
add_library (MDMTrace SHARED src/MDMTrace.cpp src/RAYTKIN1.F)
target_include_directories(MDMTrace PUBLIC ${PROJECT_SOURCE_DIR}/include)

# build executable file ${PROJECT_NAME} (MDMTracking)
add_executable(${PROJECT_NAME} MDMTraceExample.C)
target_include_directories(${PROJECT_NAME} PRIVATE  ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} PRIVATE MDMTrace JsonCpp)

# copy rayin.dat from dat/ to build/
set(DAT_FILES rayin.dat)
foreach(_datFiles ${DAT_FILES})
    configure_file(${PROJECT_SOURCE_DIR}/dat/${_datFiles} ${PROJECT_BINARY_DIR}/${_datFiles} COPYONLY)
endforeach()