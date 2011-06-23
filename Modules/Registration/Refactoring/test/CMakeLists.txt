itk_module_test()
set(ITK-RegistrationRefactoringTests
  itkNewClassExampleTest.cxx
# itkTransformVirtualDomainCalculatorTest.cxx
  itkDemonsImageToImageObjectMetricTest.cxx
  itkDeformationFieldTransformRegistrationTest.cxx
  itkArray1DToDataTest.cxx
# itkObjectToObjectOptimizerBaseTest.cxx
  itkCompositeTransformTest.cxx
  itkObjectToObjectMetricTest.cxx
  itkImageToImageObjectMetricTest.cxx
#  itkNeighborhoodNormalizedCrossCorrelationImageToImageMetricTest
#  itkGradientDescentObjectOptimizerTest.cxx
)

CreateTestDriver(ITK-RegistrationRefactoring  "${ITK-RegistrationRefactoring-Test_LIBRARIES}" "${ITK-RegistrationRefactoringTests}")

itk_add_test(NAME itkNewClassExampleTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkNewClassExampleTest)

itk_add_test(NAME itkTransformVirtualDomainCalculatorTest2D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 2 ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm)

itk_add_test(NAME itkTransformVirtualDomainCalculatorTest3D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd)

itk_add_test(NAME itkDemonsImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkDemonsImageToImageObjectMetricTest)

itk_add_test(NAME itkDemonsImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkDemonsImageToImageObjectMetricTest ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd 2)

itk_add_test(NAME itkDeformationFieldTransformRegistrationTest
          COMMAND ITK-RegistrationRefactoringTestDriver
            itkDeformationFieldTransformRegistrationTest)

itk_add_test(NAME itkArray1DToDataTest
         COMMAND ITK-RegistrationRefactoringTestDriver itkArray1DToDataTest)

itk_add_test(NAME itkObjectToObjectOptimizerBaseTest
         COMMAND ITK-RegistrationRefactoringTestDriver
         itkObjectToObjectOptimizerBaseTest)

itk_add_test(NAME itkCompositeTransformTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkCompositeTransformTest)

itk_add_test(NAME itkObjectToObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkObjectToObjectMetricTest)

itk_add_test(NAME itkImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkImageToImageObjectMetricTest)

itk_add_test(NAME itkNeighborhoodNormalizedCrossCorrelationImageToImageMetricTest)
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkNeighborhoodNormalizedCrossCorrelationImageToImageMetricTest)

itk_add_test(NAME itkGradientDescentObjectOptimizerTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkGradientDescentObjectOptimizerTest)
