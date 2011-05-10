itk_module(ITK-SpatialObjects DEPENDS  ITK-RegistrationRefactoring ITK-ImageFunction ITK-Mesh ITK-IO-Base TEST_DEPENDS ITK-TestKernel)

# Extra dependecy on RegistrationRefactoring is introduced by adding TransformParameters to TransformBase.
