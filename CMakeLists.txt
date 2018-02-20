cmake_minimum_required(VERSION 2.8)

project(Film2DoseOptimization)

find_package(OpenCV REQUIRED)

# Requires C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${OpenCV_INCLUDE_DIRS})

#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)

set(${PROJECT_SOURCE_DIR}_SCRIPTS
   README.md
   licence.txt
   step_calib_ebt3005.tif
  )

foreach(_script ${${PROJECT_SOURCE_DIR}_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


add_executable(${PROJECT_NAME} "main.cxx" ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
#target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
