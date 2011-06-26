set(DOCUMENTATION "The QuadEdgeMesh module contain a specialized set of Mesh
classes intended to represent 2-manifolds embedded in a nD space. This family
of classes provides a consistent representation of oriented surfaces and
therefore they are used as a base for implementing common mesh filters and
operations. They are commonly used for representing the output of image
segmentation algorithms.")

itk_module(ITK-QuadEdgeMesh
  DEPENDS
    ITK-Mesh
  TEST_DEPENDS
    ITK-TestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
