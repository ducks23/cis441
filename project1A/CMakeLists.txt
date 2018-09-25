cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(differencer)
SET(VTK_DIR /home/jesse/projects/VTK-build)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project1A project1A)

if(VTK_LIBRARIES)
  target_link_libraries(project1A ${VTK_LIBRARIES})
else()
  target_link_libraries(project1A vtkHybrid)
endif()

