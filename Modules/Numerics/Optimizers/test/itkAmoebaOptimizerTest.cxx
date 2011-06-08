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



#include "itkAmoebaOptimizer.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_math.h"
#include <iostream>
#include <cstdlib>



/**
 *  The objectif function is the quadratic form:
 *
 *  1/2 x^T A x - b^T x
 *
 *  Where A is represented as an itkMatrix and
 *  b is represented as a itkVector
 *
 *  The system in this example is:
 *
 *     | 3  2 ||x|   | 2|   |0|
 *     | 2  6 ||y| + |-8| = |0|
 *
 *
 *   the solution is the vector | 2 -2 |
 *
 *   and the expected final value of the function is 10.0
 *
 */
class amoebaCostFunction : public itk::SingleValuedCostFunction
{
public:

  typedef amoebaCostFunction                    Self;
  typedef itk::SingleValuedCostFunction     Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( amoebaCostFunction, SingleValuedCostFunction );

  enum { SpaceDimension=2 };

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType;

  typedef vnl_vector<double>                      VectorType;
  typedef vnl_matrix<double>                      MatrixType;


  amoebaCostFunction():m_A(SpaceDimension,SpaceDimension),m_b(SpaceDimension)
   {
    m_A[0][0] =  3;
    m_A[0][1] =  2;
    m_A[1][0] =  2;
    m_A[1][1] =  6;

    m_b[0]    =  2;
    m_b[1]    = -8;
    m_Negate = false;
    }

  double GetValue( const ParametersType & parameters ) const
    {

    VectorType v( parameters.Size() );
    for(unsigned int i=0; i<SpaceDimension; i++)
      {
      v[i] = parameters[i];
      }
    VectorType Av = m_A * v;
    double val = ( inner_product<double>( Av , v ) )/2.0;
    val -= inner_product< double >( m_b , v );
    if( m_Negate )
      {
      val *= -1.0;
      }
    return val;
    }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType & derivative ) const
    {

    VectorType v( parameters.Size() );
    for(unsigned int i=0; i<SpaceDimension; i++)
      {
      v[i] = parameters[i];
      }
    std::cout << "GetDerivative( " << v << " ) = ";
    VectorType gradient = m_A * v  - m_b;
    std::cout << gradient << std::endl;
    derivative = DerivativeType(SpaceDimension);
    for(unsigned int i=0; i<SpaceDimension; i++)
      {
      if( !m_Negate )
        {
        derivative[i] = gradient[i];
        }
      else
        {
        derivative[i] = -gradient[i];
        }
      }
    }

  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }

  // Used to switch between maximization and minimization.
  void SetNegate(bool flag )
    {
    m_Negate = flag;
    }

private:
  MatrixType        m_A;
  VectorType        m_b;
  bool              m_Negate;
};

class CommandIterationUpdateAmoeba : public itk::Command
{
public:
  typedef  CommandIterationUpdateAmoeba   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdateAmoeba()
  {
    m_IterationNumber=0;
  }
public:
  typedef itk::AmoebaOptimizer         OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( m_FunctionEvent.CheckEvent( &event ) )
        {
        std::cout << m_IterationNumber++ << "   ";
        std::cout << optimizer->GetCachedValue() << "   ";
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
        std::cout << optimizer->GetStopConditionDescription() << std::endl;
        }
      else if( m_GradientEvent.CheckEvent( &event ) )
        {
        std::cout << "Gradient " << optimizer->GetCachedDerivative() << "   ";
        }

    }
private:
  unsigned long m_IterationNumber;

  itk::FunctionEvaluationIterationEvent m_FunctionEvent;
  itk::GradientEvaluationIterationEvent m_GradientEvent;
};



int itkAmoebaOptimizerTest(int, char* [] )
{

  std::cout << "Amoeba Optimizer Test \n \n";

  typedef  itk::AmoebaOptimizer  OptimizerType;

  typedef  OptimizerType::InternalOptimizerType  vnlOptimizerType;

  // Declaration of a itkOptimizer
  OptimizerType::Pointer  itkOptimizer = OptimizerType::New();

  // set optimizer parameters
  itkOptimizer->SetMaximumNumberOfIterations( 10 );

  double xTolerance = 0.01;
  itkOptimizer->SetParametersConvergenceTolerance( xTolerance );

  double fTolerance = 0.001;
  itkOptimizer->SetFunctionConvergenceTolerance( fTolerance );

  // Declaration of the CostFunction adaptor
  amoebaCostFunction::Pointer costFunction = amoebaCostFunction::New();


  itkOptimizer->SetCostFunction( costFunction.GetPointer() );
  std::cout << "itkOptimizer->GetCostFunction(): " << itkOptimizer->GetCostFunction() << std::endl;

  vnlOptimizerType * vnlOptimizer = itkOptimizer->GetOptimizer();


  OptimizerType::ParametersType initialValue(2);       // constructor requires vector size

  initialValue[0] =  100;             // We start not far from  | 2 -2 |
  initialValue[1] = -100;

  OptimizerType::ParametersType currentValue(2);

  currentValue = initialValue;

  itkOptimizer->SetInitialPosition( currentValue );


  try
    {

    vnlOptimizer->verbose = true;

    std::cout << "Run for " << itkOptimizer->GetMaximumNumberOfIterations();
    std::cout << " iterations." << std::endl;

    itkOptimizer->StartOptimization();


    std::cout << "Continue for " << itkOptimizer->GetMaximumNumberOfIterations();
    std::cout << " iterations." << std::endl;

    itkOptimizer->SetMaximumNumberOfIterations( 100 );
    itkOptimizer->SetInitialPosition( itkOptimizer->GetCurrentPosition() );
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


  std::cout << "Number of evals = " << vnlOptimizer->get_num_evaluations() << std::endl;

  std::cout << "Optimizer: " << itkOptimizer;

  //
  // check results to see if it is within range
  //

  OptimizerType::ParametersType finalPosition;
  finalPosition = itkOptimizer->GetCurrentPosition();

  double trueParameters[2] = { 2, -2 };
  bool pass = true;

  std::cout << "Right answer   = " << trueParameters[0] << " , " << trueParameters[1] << std::endl;
  std::cout << "Final position = " << finalPosition     << std::endl;

  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > xTolerance )
      pass = false;
    }

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  // Get the final value of the optimizer
  std::cout << "Testing GetValue() : ";
  OptimizerType::MeasureType finalValue = itkOptimizer->GetValue();
  if(vcl_fabs(finalValue+9.99998)>0.01)
    {
    std::cout << "[FAILURE]" << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    std::cout << "[SUCCESS]" << std::endl;
    }


  // Set now the function to maximize
  //
  { // add a block-scope to have local variables

  std::cout << "Testing Maximization " << std::endl;

  currentValue = initialValue;

  itkOptimizer->SetInitialPosition( currentValue );

  CommandIterationUpdateAmoeba::Pointer observer =
    CommandIterationUpdateAmoeba::New();
  itkOptimizer->AddObserver( itk::IterationEvent(), observer );
  itkOptimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );


  try
    {
    // These two following statement should compensate each other
    // and allow us to get to the same result as the test above.
    costFunction->SetNegate(true);
    itkOptimizer->MaximizeOn();

    vnlOptimizer->verbose = true;

    std::cout << "Run for " << itkOptimizer->GetMaximumNumberOfIterations();
    std::cout << " iterations." << std::endl;

    itkOptimizer->StartOptimization();

    std::cout << "Continue for " << itkOptimizer->GetMaximumNumberOfIterations();
    std::cout << " iterations." << std::endl;

    itkOptimizer->SetMaximumNumberOfIterations( 100 );
    itkOptimizer->SetInitialPosition( itkOptimizer->GetCurrentPosition() );
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

  finalPosition = itkOptimizer->GetCurrentPosition();
  std::cout << "Right answer   = " << trueParameters[0] << " , " << trueParameters[1] << std::endl;
  std::cout << "Final position = " << finalPosition     << std::endl;

  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > xTolerance )
      pass = false;
    }

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  // Get the final value of the optimizer
  std::cout << "Testing GetValue() : ";
  finalValue = itkOptimizer->GetValue();
  if(vcl_fabs(finalValue+9.99998)>0.01)
    {
    std::cout << "[FAILURE]" << std::endl;
    return EXIT_FAILURE;
    }
  else
    {
    std::cout << "[SUCCESS]" << std::endl;
    }
  }




  std::cout << "Test done." << std::endl;
  return EXIT_SUCCESS;
}



