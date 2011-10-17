set(DOCUMENTATION "This module contains classes for reading and writing XML
files with the expat library.")

itk_module(ITKIOXML
  DEPENDS
    ITKExpat
    ITKIOImageBase
  TEST_DEPENDS
    ITKTestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
