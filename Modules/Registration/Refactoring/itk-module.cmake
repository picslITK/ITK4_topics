itk_module(ITKRegistrationRefactoring
  DEPENDS
    ITKRegistrationCommon
    ITKOptimizers
    ITKImageIntensity
    ITKImageFunction
    ITKImageGrid
    ITKSpatialObjects
    ITKSmoothing
    ITKImageGradient
    ITKImageFeature
    ITKFiniteDifference
    ITKHighDimensionalOptimizers
    ITKHighDimensionalMetrics
    ITKDisplacementField
  TEST_DEPENDS
    ITKTestKernel
)
