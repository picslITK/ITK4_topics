set(ITK-RegistrationRefactoring_SRC
itkGradientDescentObjectOptimizerBase.cxx
itkGradientDescentObjectOptimizer.cxx
itkQuasiNewtonObjectOptimizer.cxx
)

add_library(ITK-RegistrationRefactoring ${ITK-RegistrationRefactoring_SRC})
target_link_libraries(ITK-RegistrationRefactoring
  ${ITK-RegistrationCommon_LIBRARIES}
  ${ITK-Optimizers_LIBRARIES}
  ${ITK-ImageIntensity_LIBRARIES}
  ${ITK-ImageFunction_LIBRARIES}
  ${ITK-ImageGrid_LIBRARIES}
  ${ITK-SpatialObjects_LIBRARIES}
  ${ITK-Smoothing_LIBRARIES}
  ${ITK-ImageGradient_LIBRARIES}
  ${ITK-ImageFeature_LIBRARIES}
  ${ITK-FiniteDifference_LIBRARIES}
   )
itk_module_target(ITK-RegistrationRefactoring)
