set(ITKHighDimensionalOptimizers_SRC
itkGradientDescentObjectOptimizerBase.cxx
itkGradientDescentObjectOptimizer.cxx
)

add_library(ITKHighDimensionalOptimizers ${ITKHighDimensionalOptimizers_SRC})
target_link_libraries(ITKHighDimensionalOptimizers
  ${ITKRegistrationCommon_LIBRARIES}
  ${ITKOptimizers_LIBRARIES}
  )
itk_module_target(ITKHighDimensionalOptimizers)
