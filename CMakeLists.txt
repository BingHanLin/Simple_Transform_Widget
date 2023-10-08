cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(main)

find_package(Qt5 COMPONENTS Core Gui REQUIRED)
find_package(VTK REQUIRED)

add_executable(${PROJECT_NAME}  
  ${PROJECT_SOURCE_DIR}/main.cpp
  ${PROJECT_SOURCE_DIR}/transformRepresentation.cpp
  ${PROJECT_SOURCE_DIR}/transformWidget.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}/
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
  ${VTK_LIBRARIES}
)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS ${PROJECT_NAME}
  MODULES ${VTK_LIBRARIES}
)
