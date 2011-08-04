itk_module_test()
set(ITKRegistrationRefactoringTests
# itkTransformVirtualDomainCalculatorTest.cxx
  itkDemonsImageToImageObjectMetricTest.cxx
  itkDemonsVectorImageToVectorImageObjectMetricTest.cxx
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
  itkPolyAffineTransformTest.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2.cxx
)

set(INPUTDATA ${ITK_DATA_ROOT}/Input)
set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(ITKRegistrationRefactoring  "${ITKRegistrationRefactoring-Test_LIBRARIES}" "${ITKRegistrationRefactoringTests}")

#itk_add_test(NAME itkTransformVirtualDomainCalculatorTest2D
#      COMMAND ITKRegistrationRefactoringTestDriver #itkTransformVirtualDomainCalculatorTest 2 ${ITK_DATA_ROOT}/Input/#itkGDCMImageIOTest3.dcm ${ITK_DATA_ROOT}/Input/itkGDCMImageIOTest3.dcm)

#itk_add_test(NAME itkTransformVirtualDomainCalculatorTest3D
#      COMMAND ITKRegistrationRefactoringTestDriver #itkTransformVirtualDomainCalculatorTest 3 ${ITK_DATA_ROOT}/Input/#HeadMRVolumeWithDirection002.mhd ${ITK_DATA_ROOT}/Input/#HeadMRVolumeWithDirection002.mhd)

itk_add_test(NAME itkDemonsImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkDemonsImageToImageObjectMetricTest)

itk_add_test(NAME itkDemonsVectorImageToVectorImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkDemonsVectorImageToVectorImageObjectMetricTest)

itk_add_test(NAME itkArray1DToDataTest
         COMMAND ITKRegistrationRefactoringTestDriver itkArray1DToDataTest)

itk_add_test(NAME itkObjectToObjectOptimizerBaseTest
         COMMAND ITKRegistrationRefactoringTestDriver
         itkObjectToObjectOptimizerBaseTest)

itk_add_test(NAME itkObjectToObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver itkObjectToObjectMetricTest)

itk_add_test(NAME itkImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkImageToImageObjectMetricTest)

itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest)

itk_add_test(NAME itkGradientDescentObjectOptimizerTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkGradientDescentObjectOptimizerTest)

# NOTE - this is waiting for new testing data images to be added
#itk_add_test(NAME itkDemonsImageToImageObjectRegistrationTest
#      COMMAND ITKRegistrationRefactoringTestDriver
#      --compare ${TEMP}/itkDemonsImageToImageObjectRegistrationTest_out.nii.gz
#                DATA{Baseline/#itkDemonsImageToImageObjectRegistrationTest_r62r85_Out.nii.gz}
#      itkDemonsImageToImageObjectRegistrationTest
#      DATA{Input/r62slice.nii} DATA{Input/r85slice.nii}
#      ${TEMP}/itkDemonsImageToImageObjectRegistrationTest_out.nii.gz)

#itk_add_test(NAME itkDemonsImageToImageObjectRegistrationTest2
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkDemonsImageToImageObjectRegistrationTest2)

 itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest)

itk_add_test(NAME itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest)

itk_add_test(NAME itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest)

itk_add_test(NAME itkQuasiNewtonDemonsRegistrationTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkQuasiNewtonDemonsRegistrationTest
      ${ITK_SOURCE_DIR}/Examples/Data/r62slice.nii
      ${ITK_SOURCE_DIR}/Examples/Data/r85slice.nii
      ${TEMP}/QNDemonsOutput.nii
      )

#itk_add_test(NAME itkQuasiNewtonObjectOptimizerTest
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkQuasiNewtonObjectOptimizerTest
#      ${ITK_SOURCE_DIR}/Examples/Data/r62slice.nii
#      ${ITK_SOURCE_DIR}/Examples/Data/r85slice.nii
#      ${TEMP}/QNOptOutput.nii
#      )

itk_add_test(NAME itkPolyAffineTransformTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkPolyAffineTransformTest
      )

#TODO
# itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2)



#itk_add_test(NAME itkEuclideanDistancePointSetMetricTest
#      COMMAND ITKRegistrationRefactoringTestDriver itkPointSetToPointSetMetricTest 1)

#itk_add_test(NAME itkExpectationBasedPointSetMetricTest
#      COMMAND ITKRegistrationRefactoringTestDriver itkPointSetToPointSetMetricTest 2)

#itk_add_test(NAME itkJensenHavrdaCharvatTsallisPointSetMetricTest
#      COMMAND ITKRegistrationRefactoringTestDriver itkPointSetToPointSetMetricTest 3)
