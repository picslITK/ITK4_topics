itk_module_test()
set(ITKRegistrationRefactoringTests
  itkANTSSparseNeighborhoodCorrelationImageToImageObjectRegistrationTest.cxx
  itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest.cxx
  itkQuasiNewtonRegistrationTest.cxx
  itkQuasiNewtonDemonsRegistrationTest.cxx
  itkQuasiNewtonObjectOptimizerTest.cxx
  itkPolyAffineTransformTest.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest3.cxx
  itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest4.cxx
  # On hold until class is updated
  #itkDemonsVectorImageToVectorImageObjectRegistrationTest.cxx
  #itkDemonsVectorImageToVectorImageObjectMetricTest.cxx
)

set(INPUTDATA ${ITK_DATA_ROOT}/Input)
set(BASELINE_ROOT ${ITK_DATA_ROOT}/Baseline)
set(TEMP ${ITK_TEST_OUTPUT_DIR})

CreateTestDriver(ITKRegistrationRefactoring  "${ITKRegistrationRefactoring-Test_LIBRARIES}" "${ITKRegistrationRefactoringTests}")

# On hold until class is updated
#itk_add_test(NAME itkDemonsVectorImageToVectorImageObjectMetricTest
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkDemonsVectorImageToVectorImageObjectMetricTest)

itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkANTSNeighborhoodCorrelationImageToImageObjectMetricTest)

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

# itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest3
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest3)


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

itk_add_test(NAME itkQuasiNewtonRegistrationTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkQuasiNewtonRegistrationTest
      #${ITK_SOURCE_DIR}/Examples/Data/r62slice.nii
      #${ITK_SOURCE_DIR}/Examples/Data/r85slice.nii
      #${TEMP}/QNOptOutput.nii
      )

itk_add_test(NAME itkPolyAffineTransformTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkPolyAffineTransformTest
      )

itk_add_test(NAME itkQuasiNewtonObjectOptimizerTest
      COMMAND ITKRegistrationRefactoringTestDriver
      itkQuasiNewtonObjectOptimizerTest)

#TODO
# itk_add_test(NAME itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2
#      COMMAND ITKRegistrationRefactoringTestDriver
#      itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2)
