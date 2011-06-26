set(DOCUMENTATION "ImageFunctions are typically used as helper classes in the
implementation of ITK filters. In many cases they may act as pluggable
components at run time, while in other cases they may be implementing the
actions that a filter will apply at every iteration, or the operation a
filter may apply to every pixel. These classes are rarely used from
applications. Their main role is to assist in the implementation of ITK
filters.")

itk_module(ITK-ImageFunction
  DEPENDS
    ITK-ImageAdaptors
    ITK-Transform
    ITK-Statistics
  TEST_DEPENDS
    ITK-TestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
