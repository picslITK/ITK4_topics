
cmake_minimum_required(VERSION 2.4)
#Change PROJECT_NAME to the name of your project
PROJECT(identity_test)

# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 


SET(IDENTITY_TEST_INCLUDE_DIRS
/Users/brianavants/code/eigen/Eigen
)
SET(IDENTITY_TEST_SRC_DIR
./
)

INCLUDE_DIRECTORIES(${IDENTITY_TEST_INCLUDE_DIRS})
LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )


ADD_EXECUTABLE(identity_test ${IDENTITY_TEST_SRC_DIR}identity_test.cxx)
TARGET_LINK_LIBRARIES(identity_test ITKCommon  ITKIO ITKNumerics  ITKStatistics )

