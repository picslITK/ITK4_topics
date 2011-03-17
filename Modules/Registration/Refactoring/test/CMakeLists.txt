itk_module_test()
set(ITK-RegistrationRefactoringTests
  itkNewClassExampleTest.cxx
)

CreateTestDriver(ITK-RegistrationRefactoring  "${ITK-RegistrationRefactoring-Test_LIBRARIES}" "${ITK-RegistrationRefactoringTests}")


add_test(NAME itkNewClassExampleTest
      COMMAND ITK-RegistrationRefactoringTestDriver itkNewClassExampleTest)
