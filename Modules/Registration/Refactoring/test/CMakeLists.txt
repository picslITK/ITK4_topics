itk_module_test()
set(ITK-RegistrationRefactoringTests
  itkNewClassExampleTest.cxx
# itkTransformVirtualDomainCalculatorTest.cxx
  itkDemonsImageToImageObjectMetricTest.cxx
  itkDeformationFieldTransformRegistrationTest.cxx
  itkArray1DToDataTest.cxx
  itkObjectToObjectOptimizerBaseTest.cxx
  itkObjectToObjectMetricTest.cxx
#  itkPointSetToPointSetMetricTest.cxx
  itkImageToImageObjectMetricTest.cxx
  itkGradientDescentObjectOptimizerTest.cxx
  itkDemonsImageToImageObjectRegistrationTest.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest.cxx
#  itkDemonsImageToImageObjectRegistrationTest2.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest.cxx
  itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest.cxx
  itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest.cxx
  itkQuasiNewtonObjectOptimizerTest.cxx
  itkQuasiNewtonDemonsRegistrationTest.cxx
)

set(INPUTDATA ${ITK_DATA_ROOT}/Input)
set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(ITK-RegistrationRefactoring  "${ITK-RegistrationRefactoring-Test_LIBRARIES}" "${ITK-RegistrationRefactoringTests}")

itk_add_test(NAME itkNewClassExampleTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkNewClassExampleTest)

#itk_add_test(NAME itkTransformVirtualDomainCalculatorTest2D
#      COMMAND ITK-RegistrationRefactoringTestDriver #itkTransformVirtualDomainCalculatorTest 2 ${ITK_DATA_ROOT}/Input/#itkGDCMImageIOTest3.dcm ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm)

#itk_add_test(NAME itkTransformVirtualDomainCalculatorTest3D
#      COMMAND ITK-RegistrationRefactoringTestDriver #itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/#HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/#HeadMRVolumeWithDirection002.mhd)

itk_add_test(NAME itkDemonsImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkDemonsImageToImageObjectMetricTest)

itk_add_test(NAME itkDeformationFieldTransformRegistrationTest
          COMMAND ITK-RegistrationRefactoringTestDriver
            itkDeformationFieldTransformRegistrationTest)

itk_add_test(NAME itkArray1DToDataTest
         COMMAND ITK-RegistrationRefactoringTestDriver itkArray1DToDataTest)

itk_add_test(NAME itkObjectToObjectOptimizerBaseTest
         COMMAND ITK-RegistrationRefactoringTestDriver
         itkObjectToObjectOptimizerBaseTest)

itk_add_test(NAME itkObjectToObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkObjectToObjectMetricTest)

itk_add_test(NAME itkImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkImageToImageObjectMetricTest)

itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest)

itk_add_test(NAME itkGradientDescentObjectOptimizerTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkGradientDescentObjectOptimizerTest)

itk_add_test(NAME itkDemonsImageToImageObjectRegistrationTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      --compare ${TEMP}/itkDemonsImageToImageObjectRegistrationTest_out.nii.gz
                DATA{Baseline/itkDemonsImageToImageObjectRegistrationTest_r62r85_Out.nii.gz}
      itkDemonsImageToImageObjectRegistrationTest
      DATA{Input/r62slice.nii} DATA{Input/r85slice.nii}
      ${TEMP}/itkDemonsImageToImageObjectRegistrationTest_out.nii.gz)

#itk_add_test(NAME itkDemonsImageToImageObjectRegistrationTest2
#      COMMAND ITK-RegistrationRefactoringTestDriver
#      itkDemonsImageToImageObjectRegistrationTest2)

itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest)

itk_add_test(NAME itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest)

itk_add_test(NAME itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest)

itk_add_test(NAME itkQuasiNewtonDemonsRegistrationTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkQuasiNewtonDemonsRegistrationTest
      ${ITK_SOURCE_DIR}/Examples/Data/r62slice.nii
      ${ITK_SOURCE_DIR}/Examples/Data/r85slice.nii
      ${TEMP}/QNDemonsOutput.nii
      )

itk_add_test(NAME itkQuasiNewtonObjectOptimizerTest
      COMMAND ITK-RegistrationRefactoringTestDriver
      itkQuasiNewtonObjectOptimizerTest
      ${ITK_SOURCE_DIR}/Examples/Data/r62slice.nii
      ${ITK_SOURCE_DIR}/Examples/Data/r85slice.nii
      ${TEMP}/QNOptOutput.nii
      )

#itk_add_test(NAME itkEuclideanDistancePointSetMetricTest
#      COMMAND ITK-ReviewTestDriver itkPointSetToPointSetMetricTest 1)

#itk_add_test(NAME itkExpectationBasedPointSetMetricTest
#      COMMAND ITK-ReviewTestDriver itkPointSetToPointSetMetricTest 2)

#itk_add_test(NAME itkJensenHavrdaCharvatTsallisPointSetMetricTest
#      COMMAND ITK-ReviewTestDriver itkPointSetToPointSetMetricTest 3)
