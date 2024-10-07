cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

#set(CMAKE_GENERATOR_TOOLSET "ClangCL")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(DALPHA7X_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# add required
set(PROJECT_REQUIRED_DLL embree4.dll tbb12.dll tbbmalloc.dll)
set(PROJECT_LIBRARY_DLL_DIR ${DALPHA7X_ROOT_DIR}/thirdparty/embree/dll)

project(AlphaX7Renderer VERSION 0.1.0)

# source file
file(GLOB_RECURSE HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.h")
file(GLOB_RECURSE HPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.hpp")
file(GLOB_RECURSE SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")
file(GLOB_RECURSE SOURCE_C_FILES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.c")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/Source" FILES ${HEADER_FILES} ${SOURCE_FILES} ${SOURCE_C_FILES} ${HPP_FILES})

# include source directories
include_directories(${DALPHA7X_ROOT_DIR}/source)

# add thirdparty
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)

set (ALL_PBRT_LIBS
  double-conversion
)

# link thirdparty libs
# link double conversion
link_libraries (${ALL_PBRT_LIBS})
# link embree
link_libraries(${DALPHA7X_ROOT_DIR}/thirdparty/embree/lib/embree4.lib)
link_libraries(${DALPHA7X_ROOT_DIR}/thirdparty/embree/lib/tbb12.lib)

# include thiryparty
# include embree headers
include_directories(${DALPHA7X_ROOT_DIR}/thirdparty/embree/include)
# include double conversion
include_directories(${DOUBLE_CONVERSION_INCLUDE})
# include gmtl
include_directories(${GLM_INCLUDE})



# defines
add_definitions(-DALPHA7X_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# executable
add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES} ${SOURCE_C_FILES})

# copy embree dll
foreach (REQUIRED_DLL  ${PROJECT_REQUIRED_DLL})
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${PROJECT_LIBRARY_DLL_DIR}/${REQUIRED_DLL}"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>")
        message(STATUS "COPY ST_DLL ${REQUIRED_DLL}")
endforeach (REQUIRED_DLL)
