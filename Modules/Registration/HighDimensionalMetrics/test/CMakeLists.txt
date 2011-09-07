itk_module_test()
set(ITKHighDimensionalMetricsTests
  itkObjectToObjectMetricTest.cxx
  itkImageToImageObjectMetricTest.cxx
  itkDemonsImageToImageObjectMetricTest.cxx
  itkMattesMutualInformationImageToImageObjectMetricTest.cxx
)

set(INPUTDATA ${ITK_DATA_ROOT}/Input)
set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(ITKHighDimensionalMetrics  "${ITKHighDimensionalMetrics-Test_LIBRARIES}" "${ITKHighDimensionalMetricsTests}")

itk_add_test(NAME itkObjectToObjectMetricTest
      COMMAND ITKHighDimensionalMetricsTestDriver itkObjectToObjectMetricTest)

itk_add_test(NAME itkImageToImageObjectMetricTest
      COMMAND ITKHighDimensionalMetricsTestDriver
              itkImageToImageObjectMetricTest)

itk_add_test(NAME itkDemonsImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkDemonsImageToImageObjectMetricTest)

itk_add_test(NAME itkMattesMutualInformationImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkMattesMutualInformationImageToImageObjectMetricTest)
