/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkLBFGSBOptimizer.h"
#include "itkTextOutput.h"
#include "vnl/vnl_math.h"
#include <iostream>

/**
 *  The objective function is the quadratic form:
 *
 *  f(x) = 1/2 x^T A x - b^T x  subject to  -1 <= x <= 10
 *
 *  Where A is represented as an itkMatrix and
 *  b is represented as a itkVector
 *
 *  The system in this example is:
 *
 *          | 3  2 |       | 2|
 *      A=  | 2  6 |   b = |-8|
 *
 *   the solution is the vector | 4/3 -1 |
 *
 */
class LBFGSBCostFunction : public itk::SingleValuedCostFunction
{
public:

  typedef LBFGSBCostFunction                    Self;
  typedef itk::SingleValuedCostFunction     Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( LBFGSBCostFunction, SingleValuedCostFunction );

  enum { SpaceDimension=2 };

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;

  typedef vnl_vector<double>                      VectorType;
  typedef vnl_matrix<double>                      MatrixType;

  typedef double MeasureType ;

  LBFGSBCostFunction()
  {
  }

  double GetValue( const ParametersType & position ) const
  {

    double x = position[0];
    double y = position[1];

    std::cout << "GetValue ( " ;
    std::cout << x << " , " << y;
    std::cout << ") = ";

    double val = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y;

    std::cout << val << std::endl;

    return val;
  }

  void GetDerivative( const ParametersType & position,
                            DerivativeType  & derivative ) const
  {

    double x = position[0];
    double y = position[1];

    std::cout << "GetDerivative ( " ;
    std::cout << x << " , " << y;
    std::cout << ") = ";

    derivative = DerivativeType(SpaceDimension);
    derivative[0] = 3*x + 2*y -2;
    derivative[1] = 2*x + 6*y +8;
    std::cout << "(" ;
    std::cout << derivative[0] <<" , ";
    std::cout << derivative[1] << ")" << std::endl;

  }


  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }

private:


};



int itkLBFGSBOptimizerTest(int, char *[])
{

  itk::OutputWindow::SetInstance(itk::TextOutput::New().GetPointer());

  std::cout << "LBFGSB Optimizer Test \n \n";

  typedef  itk::LBFGSBOptimizer  OptimizerType;

  // Declaration of a itkOptimizer
  OptimizerType::Pointer  itkOptimizer = OptimizerType::New();

  itkOptimizer->Print( std::cout );


  // Declaration of the CostFunction adaptor
  LBFGSBCostFunction::Pointer costFunction = LBFGSBCostFunction::New();


  itkOptimizer->SetCostFunction( costFunction.GetPointer() );

  const double F_Convergence_Factor  = 1e+7;      // Function value tolerance
  const double Projected_G_Tolerance = 1e-5;      // Proj gradient tolerance
  const int    Max_Iterations   =   100; // Maximum number of iterations

  itkOptimizer->SetCostFunctionConvergenceFactor( F_Convergence_Factor );
  itkOptimizer->SetProjectedGradientTolerance( Projected_G_Tolerance );
  itkOptimizer->SetMaximumNumberOfIterations( Max_Iterations );
  itkOptimizer->SetMaximumNumberOfEvaluations( Max_Iterations );

  const unsigned int SpaceDimension = 2;
  OptimizerType::ParametersType initialValue(SpaceDimension);

  // Starting point
  initialValue[0] =  10;
  initialValue[1] =  10;

  OptimizerType::ParametersType currentValue(2);

  currentValue = initialValue;

  itkOptimizer->SetInitialPosition( currentValue );

  // Set up boundary conditions
  OptimizerType::BoundValueType lower(SpaceDimension);
  OptimizerType::BoundValueType upper(SpaceDimension);
  OptimizerType::BoundSelectionType select(SpaceDimension);

  lower.Fill( -1 );
  upper.Fill( 10 );
  select.Fill( 2 );

  itkOptimizer->SetLowerBound( lower );
  itkOptimizer->SetUpperBound( upper );
  itkOptimizer->SetBoundSelection( select );

  itkOptimizer->Print( std::cout );

  try
    {

    itkOptimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  const OptimizerType::ParametersType & finalPosition = itkOptimizer->GetCurrentPosition();

  std::cout << "Solution        = ("
    << finalPosition[0] << ","
    << finalPosition[1] << ")" << std::endl;
  std::cout << "Final Function Value = "
    << itkOptimizer->GetValue() << std::endl;

  std::cout << "Infinity Norm of Projected Gradient = "
    << itkOptimizer->GetInfinityNormOfProjectedGradient() << std::endl;
  std::cout << "End condition   = "
    << itkOptimizer->GetStopConditionDescription() << std::endl;
  std::cout << "Trace   = " << itkOptimizer->GetTrace() << std::endl;
  std::cout << "CostFunctionConvergenceFactor   = "
    << itkOptimizer->GetCostFunctionConvergenceFactor() << std::endl;
  std::cout << "ProjectedGradientTolerance   = "
    << itkOptimizer->GetProjectedGradientTolerance() << std::endl;
  std::cout << "MaximumNumberOfIterations   = "
    << itkOptimizer->GetMaximumNumberOfIterations() << std::endl;
  std::cout << "MaximumNumberOfEvaluations   = "
    << itkOptimizer->GetMaximumNumberOfEvaluations() << std::endl;
  std::cout << "MaximumNumberOfCorrections   = "
    << itkOptimizer->GetMaximumNumberOfCorrections() << std::endl;

  //
  // check results to see if it is within range
  //
  bool pass = true;
  std::string errorIn;

  double trueParameters[2] = { 4.0/3.0, -1.0 };
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      {
      pass = false;
      errorIn = "solution";
      }
    }

  if( vnl_math_abs( itkOptimizer->GetValue() - -7.66667 ) > 0.01 )
    {
    pass = false;
    errorIn = "final function value";
    }


  if( vnl_math_abs( itkOptimizer->GetInfinityNormOfProjectedGradient()
      -  1.77636e-15 ) > 0.01 )
    {
    pass = false;
    errorIn = "infinity norm of projected gradient";
    }

  if( !pass )
    {
    std::cout << "\nError in " << errorIn << ".\n";
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;

}
