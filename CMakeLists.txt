cmake_minimum_required(VERSION 2.8)
project(fluid-simulation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(GLFW_SOURCE_DIR "external/glfw")
set(GLM_SOURCE_DIR "external/glm")
set(GLEW_SOURCE_DIR "external/glew")
set(LODEPNG_SOURCE_DIR "external/lodepng")

find_package(OpenGL REQUIRED)


# Building only the GLFW lib
set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "")
set(GLFW_BUILD_TESTS OFF CACHE BOOL "")
set(GLFW_BUILD_DOCS OFF CACHE BOOL "")
set(GLFW_INSTALL OFF CACHE BOOL "")

add_subdirectory(${GLFW_SOURCE_DIR} ${GLEW_SOURCE_DIR})

include_directories(${GLFW_SOURCE_DIR}/include
        ${GLFW_SOURCE_DIR}/deps
        ${GLM_SOURCE_DIR}
        ${GLEW_SOURCE_DIR}/include)
include_directories(${OPENGL_INCLUDE_DIRS})

add_definitions(
        -DGLEW_STATIC
)

add_library(glew STATIC
        ${GLEW_SOURCE_DIR}/src/glew.c
        ${GLEW_SOURCE_DIR}/include
        )

ADD_LIBRARY(lodepng
        ${LODEPNG_SOURCE_DIR}/lodepng.cpp
        ${LODEPNG_SOURCE_DIR}/lodepng.h
        )


target_link_libraries(glew ${OPENGL_LIBRARIES})



add_executable(fluid-simulation main.cpp Renderer.cpp Renderer.h FluidSimulator.cpp FluidSimulator.h Shader.cpp Shader.h)
target_link_libraries("fluid-simulation" lodepng glfw ${OPENGL_LIBRARIES} ${GLFW_LIBRARIES} glew)

#if (APPLE)
#  set(APP_NAME "OpenGLBoilerplate")
#
#  if (EXISTS icon.icns)
#    add_executable(${APP_NAME} MACOSX_BUNDLE main.cpp icon.icns)
#    set_source_files_properties("icon.icns" PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
#  else ()
#    add_executable(${APP_NAME} MACOSX_BUNDLE main.cpp)
#  endif ()
#
#  target_link_libraries(${APP_NAME} glfw ${OPENGL_glu_LIBRARY} ${GLFW_LIBRARIES} glew)
#  set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_ICON_FILE "icon.icns")
#
#  set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME ${APP_NAME})
#  set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_SHORT_VERSION_STRING "1.0")
#  set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_LONG_VERSION_STRING "1.0")
#  set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_INFO_STRING "OpenGL boilerplate example app")
#endif()