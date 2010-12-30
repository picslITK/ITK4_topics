
cmake_minimum_required(VERSION 2.4)
#Change PROJECT_NAME to the name of your project
PROJECT(BostonCodeParty)

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


SET(IDENTITY_TEST_SRC_DIR
./
)

INCLUDE_DIRECTORIES(${IDENTITY_TEST_INCLUDE_DIRS})
LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )

ADD_EXECUTABLE(sym_metric ${IDENTITY_TEST_SRC_DIR}sym_metric.cxx)
TARGET_LINK_LIBRARIES(sym_metric ITKCommon  ITKIO ITKNumerics  ITKStatistics )

ADD_EXECUTABLE(demons ${IDENTITY_TEST_SRC_DIR}demons.cxx)
TARGET_LINK_LIBRARIES(demons ITKCommon  ITKIO ITKNumerics  ITKStatistics )

ADD_EXECUTABLE(affine ${IDENTITY_TEST_SRC_DIR}affine.cxx)
TARGET_LINK_LIBRARIES(affine ITKCommon  ITKIO ITKNumerics  ITKStatistics )

ADD_EXECUTABLE(templated_resampler ${IDENTITY_TEST_SRC_DIR}templated_resampler.cxx)
TARGET_LINK_LIBRARIES(templated_resampler ITKCommon  ITKIO ITKNumerics  ITKStatistics )


# 1 sym(metric) derived from objecttoobjectmetric and used in registration -- both a dense and non-dense version of metric with maximal code re-use

# 2 affine registration via composite tx

# 3 affine+deformable via composite tx with sparse / dense metrics attached

# 4 test speed of resampleimagefilter if we template all the  parameters --- interpolator, transform, etc.   i imagine this would be notably faster, given the basic_functor_test.  consider templating interpolation in metrics?
