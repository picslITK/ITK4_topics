itk_module(ITK-Transform DEPENDS ITK-ImageStatistics ITK-RegistrationRefactoring ITK-Review TEST_DEPENDS ITK-TestKernel ITK-ImageFunction ITK-ImageGrid ITK-SpatialObjects)
# Extra test depedency on ImageFunction and ImageGrid is introduced by itkBSplineDeformableTransformTest.
# Extra test dependency on  SpatialObjects is introduced by itkCenteredVersorTransformInitializerTest.
# Extra dependecy on RegistrationRefactoring is introduced by adding TransformParameters to TransformBase.
# Extra dependecy on Review is introduced by reference to itkDeformationFieldTransform.h.
