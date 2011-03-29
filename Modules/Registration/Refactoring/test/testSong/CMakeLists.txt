project(CompositiveAffineProject)

cmake_minimum_required(VERSION 2.8.2)
# set(CMAKE_BUILD_TYPE "Debug")

#Change PROJECT_NAME to the name of your project

option(USE_ITK "Use ITK Libraries" ON)
option(USE_VTK "Use VTK Libraries" OFF)

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)


# Set up ITK
if(USE_ITK)
  find_package(ITK)
  if(ITK_FOUND)
    include(${ITK_USE_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  endif(ITK_FOUND)
endif(USE_ITK)


#The following lines are required to use Dart
#include(CTest)
# enable_testing()



LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )


add_executable(itkCompositiveAffineTest itkCompositiveAffineTest.cxx )
target_link_libraries(itkCompositiveAffineTest ${ITK_LIBRARIES})
