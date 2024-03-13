cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(main)

# Generate clangd configuration(for developing enviroment only)
if(MSVC)
    file(WRITE ${PROJECT_SOURCE_DIR}/.clangd "CompileFlags:\n")
    file(APPEND ${PROJECT_SOURCE_DIR}/.clangd "    CompilationDatabase: ${CMAKE_BINARY_DIR}\n")
    file(APPEND ${PROJECT_SOURCE_DIR}/.clangd "    Add:\n")

    set(MSVC_STANDARD_INCLUDE_PATHS "$ENV{INCLUDE}")
    foreach(ONE_PATH ${MSVC_STANDARD_INCLUDE_PATHS})
        file(APPEND ${PROJECT_SOURCE_DIR}/.clangd "        - \"-I${ONE_PATH}\"\n")
    endforeach()
endif()

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

vtk_module_autoinit(
  TARGETS ${PROJECT_NAME}
  MODULES ${VTK_LIBRARIES}
)

install(TARGETS  ${PROJECT_NAME}
        DESTINATION .
)