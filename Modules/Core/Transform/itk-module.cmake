set(DOCUMENTATION "Transforms are essential components of image registration
framework in ITK. They are typically used for representing the mapping between
the physical coordinate system of one image and the physical coordinate system
of another image. They are also commonly used in the process of resampling
images, particulaly when mapping them between coordinate systems. Transforms
are a large family in ITK and form a prolific group of classes in the
toolkit.")

itk_module(ITK-Transform DEPENDS ITK-Common ITK-Statistics ITK-HDF5 TEST_DEPENDS ITK-TestKernel ITK-ImageFunction ITK-ImageGrid ITK-SpatialObjects DESCRIPTION "${DOCUMENTATION}")

# Extra dependency on Common is introduced by itkBSplineDeformationFieldTransform
# Extra test depedency on ImageFunction and ImageGrid is introduced by itkBSplineDeformableTransformTest.
# Extra test dependency on  SpatialObjects is introduced by itkCenteredVersorTransformInitializerTest.
