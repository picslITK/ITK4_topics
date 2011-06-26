set(DOCUMENTATION "This module contains the base classes needed to implement
finite differences image filters. They include both sparse and dense finite
differences. Most of the classes in this module are abstract and therefore are
not intended to be instantiated by themselves.")

itk_module(ITK-FiniteDifference
  DEPENDS
    ITK-Common
  TEST_DEPENDS
    ITK-TestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
