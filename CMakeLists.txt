cmake_minimum_required(VERSION 3.10)
project(MDMTraceExample)

set(CMAKE_CXX_STANDARD 17)

enable_language (Fortran)
set (CMAKE_Fortran_FLAGS "-fPIC -finit-local-zero -falign-commons -fno-automatic")
set (CMAKE_CXX_FLAGS "-fPIC")

# build libraray libJsonCpp
add_library(JsonCpp SHARED src/jsoncpp.cpp)
target_include_directories(JsonCpp PUBLIC ${PROJECT_SOURCE_DIR}/include/json)

# build library libMDMFieldMapTransport
add_library (MDMFieldMapTransport SHARED
    src/MDMFieldMap.cpp
    src/MDMFieldMapTrace.cpp
)
target_include_directories(MDMFieldMapTransport PUBLIC ${PROJECT_SOURCE_DIR}/include)

# build library libMDMTraceFit
add_library (MDMTrace SHARED
    src/MDMTrace.cpp
    src/RAYTKIN1.F
    src/MDMFieldMapInterop.F90
)
target_include_directories(MDMTrace PUBLIC ${PROJECT_SOURCE_DIR}/include)

# build executable file ${PROJECT_NAME} (MDMTracking)
add_executable(${PROJECT_NAME} MDMTraceExample.C)
target_include_directories(${PROJECT_NAME} PRIVATE  ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} PRIVATE MDMTrace JsonCpp)

add_executable(MDMFieldMapGenerator src/MDMFieldMapGenerator.cpp)
target_include_directories(MDMFieldMapGenerator PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MDMFieldMapGenerator PRIVATE MDMTrace MDMFieldMapTransport JsonCpp)

add_executable(MDMFieldMapTraceExample MDMFieldMapTraceExample.cpp)
target_include_directories(MDMFieldMapTraceExample PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MDMFieldMapTraceExample PRIVATE MDMFieldMapTransport JsonCpp)

add_executable(GenerateIonOptics GenerateIonOptics.cpp)
target_include_directories(GenerateIonOptics PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_definitions(GenerateIonOptics PRIVATE MDMTRACE_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(GenerateIonOptics PRIVATE MDMFieldMapTransport JsonCpp)

find_package(ROOT QUIET COMPONENTS Core RIO Hist Graf Gpad Matrix)
if(ROOT_FOUND AND ROOT_LIBRARIES)
    set(MDMTRACE_ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS})
    set(MDMTRACE_ROOT_LIBRARIES ${ROOT_LIBRARIES})
    separate_arguments(MDMTRACE_ROOT_COMPILE_OPTIONS UNIX_COMMAND "${ROOT_CXX_FLAGS}")
else()
    find_program(ROOT_CONFIG_EXECUTABLE root-config)
    if(NOT ROOT_CONFIG_EXECUTABLE)
        message(FATAL_ERROR "ROOT is required to build Compare and GenerateIonOptics. Install ROOT or ensure root-config is on PATH.")
    endif()
    execute_process(
        COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
        RESULT_VARIABLE ROOT_CONFIG_CFLAGS_RESULT
        OUTPUT_VARIABLE ROOT_CONFIG_CFLAGS
        ERROR_VARIABLE ROOT_CONFIG_CFLAGS_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT ROOT_CONFIG_CFLAGS_RESULT EQUAL 0)
        message(FATAL_ERROR "root-config --cflags failed: ${ROOT_CONFIG_CFLAGS_ERROR}")
    endif()
    execute_process(
        COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
        RESULT_VARIABLE ROOT_CONFIG_LIBS_RESULT
        OUTPUT_VARIABLE ROOT_CONFIG_LIBS
        ERROR_VARIABLE ROOT_CONFIG_LIBS_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT ROOT_CONFIG_LIBS_RESULT EQUAL 0)
        message(FATAL_ERROR "root-config --libs failed: ${ROOT_CONFIG_LIBS_ERROR}")
    endif()
    separate_arguments(MDMTRACE_ROOT_COMPILE_OPTIONS UNIX_COMMAND "${ROOT_CONFIG_CFLAGS}")
    separate_arguments(MDMTRACE_ROOT_LIBRARIES UNIX_COMMAND "${ROOT_CONFIG_LIBS}")
endif()

add_executable(Compare Compare.cpp)
target_include_directories(Compare PRIVATE ${PROJECT_SOURCE_DIR}/include ${MDMTRACE_ROOT_INCLUDE_DIRS})
target_compile_options(Compare PRIVATE ${MDMTRACE_ROOT_COMPILE_OPTIONS})
target_compile_definitions(Compare PRIVATE MDMTRACE_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(Compare PRIVATE MDMTrace MDMFieldMapTransport JsonCpp ${MDMTRACE_ROOT_LIBRARIES})

target_include_directories(GenerateIonOptics PRIVATE ${MDMTRACE_ROOT_INCLUDE_DIRS})
target_compile_options(GenerateIonOptics PRIVATE ${MDMTRACE_ROOT_COMPILE_OPTIONS})
target_link_libraries(GenerateIonOptics PRIVATE ${MDMTRACE_ROOT_LIBRARIES})

# copy rayin.dat from dat/ to build/
set(DAT_FILES rayin.dat)
foreach(_datFiles ${DAT_FILES})
    configure_file(${PROJECT_SOURCE_DIR}/dat/${_datFiles} ${PROJECT_BINARY_DIR}/${_datFiles} COPYONLY)
endforeach()
