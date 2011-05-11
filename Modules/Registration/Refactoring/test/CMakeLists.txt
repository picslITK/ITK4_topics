itk_module_test()
set(ITK-RegistrationRefactoringTests
  itkNewClassExampleTest.cxx
  itkTransformVirtualDomainCalculatorTest.cxx
  itkDemonsImageToImageMetricTest.cxx
  itkDeformationFieldTransformRegistrationTest.cxx
  itkTransformParametersTest.cxx
  itkImageVectorTransformParametersHelperTest.cxx
)

CreateTestDriver(ITK-RegistrationRefactoring  "${ITK-RegistrationRefactoring-Test_LIBRARIES}" "${ITK-RegistrationRefactoringTests}")

add_test(NAME itkNewClassExampleTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkNewClassExampleTest)
add_test(NAME itkTransformVirtualDomainCalculatorTest2D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 2 ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm)
add_test(NAME itkTransformVirtualDomainCalculatorTest3D
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd)
add_test(NAME itkDemonsImageToImageMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/HeadMRVolumeWithDirection002.mhd)
add_test(NAME itkDeformationFieldTransformRegistrationTest
          COMMAND ITK-RegistrationRefactoringTestDriver
            itkDeformationFieldTransformRegistrationTest)

add_test(NAME itkTransformParametersTest
         COMMAND ITK-RegistrationRefactoringTestDriver
          itkTransformParametersTest)
add_test(NAME itkImageVectorTransformParametersHelperTest
         COMMAND ITK-RegistrationRefactoringTestDriver
          itkImageVectorTransformParametersHelperTest)
