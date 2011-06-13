itk_module_test()
set(ITK-RegistrationRefactoringTests
  itkNewClassExampleTest.cxx
#  itkTransformVirtualDomainCalculatorTest.cxx
#  itkDemonsImageToImageMetricTest.cxx
  itkDeformationFieldTransformRegistrationTest.cxx
  itkArray1DToDataTest.cxx
#  itkObjectToObjectThreadedMetricOptimizerBaseTest.cxx
  itkCompositeTransformTest.cxx
  itkObjectToObjectMetricTest.cxx
)

CreateTestDriver(ITK-RegistrationRefactoring  "${ITK-RegistrationRefactoring-Test_LIBRARIES}" "${ITK-RegistrationRefactoringTests}")

itk_add_test(NAME itkNewClassExampleTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkNewClassExampleTest)

itk_add_test(NAME itkTransformVirtualDomainCalculatorTest2D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 2 ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm)

itk_add_test(NAME itkTransformVirtualDomainCalculatorTest3D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd)

#itk_add_test(NAME itkDemonsImageToImageMetricTest
#      COMMAND ITK-RegistrationRefactoringTestDriver itkDemonsImageToImageMetricTest ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd 2)

itk_add_test(NAME itkDeformationFieldTransformRegistrationTest
          COMMAND ITK-RegistrationRefactoringTestDriver
            itkDeformationFieldTransformRegistrationTest)

itk_add_test(NAME itkArray1DToDataTest
         COMMAND ITK-RegistrationRefactoringTestDriver itkArray1DToDataTest)

itk_add_test(NAME itkObjectToObjectThreadedMetricOptimizerBaseTest
         COMMAND ITK-RegistrationRefactoringTestDriver
         itkObjectToObjectThreadedMetricOptimizerBaseTest)
itk_add_test(NAME itkCompositeTransformTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkCompositeTransformTest)
itk_add_test(NAME itkObjectToObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkObjectToObjectMetricTest)
