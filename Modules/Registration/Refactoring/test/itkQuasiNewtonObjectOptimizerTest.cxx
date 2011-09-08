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

#include "itkObjectToObjectMetric.h"
#include "itkQuasiNewtonObjectOptimizer.h"
#include "vnl/vnl_math.h"

//We need this as long as we have to define ImageToData as a fwd-declare
// in itkImageToImageObjectMetric.h
#include "itkImageToData.h"

using namespace itk;

namespace {
/**
 *  \class RosenbrockMetric for test
 *
 *  The objective function is Rosenbrock's Function, which has a banana shape.
 *
 *  f(x) = 100*(y-x*x)^2+(1-x)^2
 *
 *  The Quasi-Newton method for this function shape is much faster than a standard
 *  Gradient Descent method.
 *
 */
class RosenbrockMetric : public itk::ObjectToObjectMetric
{
public:

  typedef RosenbrockMetric                Self;
  typedef itk::ObjectToObjectMetric       Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( RosenbrockMetric, ObjectToObjectMetric );

  enum { SpaceDimension=2 };

  typedef Superclass::ParametersType        ParametersType;
  typedef Superclass::ParametersValueType   ParametersValueType;
  typedef Superclass::DerivativeType        DerivativeType;
  typedef Superclass::MeasureType           MeasureType;

  RosenbrockMetric()
  {
    m_Parameters.SetSize( SpaceDimension );
    m_Parameters.Fill( 0 );
  }

  void Initialize(void) throw ( itk::ExceptionObject ) {}

  void GetValueAndDerivative( MeasureType & value,
                              DerivativeType & derivative )
  {
    this->SetDebug(true);

    if( derivative.Size() != 2 )
      derivative.SetSize(2);

    double x = m_Parameters[0];
    double y = m_Parameters[1];

    double u = y-x*x;
    double v = 1-x;
    value = 100*u*u+v*v;
    value = -value; //for maximization

    /* The optimizer simply takes the derivative from the metric
     * and adds it to the transform after scaling. So instead of
     * setting a 'minimize' option in the gradient, we return
     * a minimizing derivative. */
    derivative[0] = 200*u*(-2*x) + 2*v*(-1);
    derivative[1] = 200*u;
    derivative[0] = -derivative[0]; //for maximization
    derivative[1] = -derivative[1];

    if (this->GetDebug())
      {
      std::cout << "GetValueAndDerivative(";
      std::cout << x << " ";
      std::cout << y << ") returns " << std::endl;
      std::cout << "  value: " << value << std::endl;
      std::cout << "  derivative: " << derivative << std::endl;
      }
  }

  MeasureType  GetValue()
  {
    return 0.0;
  }

  void UpdateTransformParameters( DerivativeType & update, ParametersValueType )
  {
    m_Parameters += update;
  }

  unsigned int GetNumberOfParameters(void) const
  {
    return SpaceDimension;
  }

  unsigned int GetNumberOfLocalParameters() const
  {
    return SpaceDimension;
  }

  bool HasLocalSupport() const
  {
    return false;
  }

  /* These Set/Get methods are only needed for this test derivation that
   * isn't using a transform */
  void SetParameters( ParametersType & parameters )
  {
    m_Parameters = parameters;
  }

  const ParametersType & GetParameters() const
  {
    return m_Parameters;
  }

private:

  ParametersType m_Parameters;
};

}//namespace

int itkQuasiNewtonObjectOptimizerTest(int, char* [] )
{
  std::cout << "QuasiNewtonObjectOptimizer Test ";
  std::cout << std::endl << std::endl;

  typedef  itk::QuasiNewtonObjectOptimizer      OptimizerType;

  typedef OptimizerType::ScalesType             ScalesType;

  // Declaration of a itkOptimizer
  OptimizerType::Pointer  itkOptimizer = OptimizerType::New();

  // Declaration of the Metric
  RosenbrockMetric::Pointer metric = RosenbrockMetric::New();

  itkOptimizer->SetMetric( metric );

  typedef RosenbrockMetric::ParametersType    ParametersType;

  const unsigned int spaceDimension =
                      metric->GetNumberOfParameters();

  ParametersType  initialPosition( spaceDimension );

  initialPosition[0] = 10;//-200;
  initialPosition[1] = -50;
  metric->SetParameters( initialPosition );

  //itkOptimizer->SetLearningRate( 0.1 );
  itkOptimizer->SetNumberOfIterations( 1000 );
  itkOptimizer->SetLineSearchEnabled(true);

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

  ParametersType finalPosition = metric->GetParameters();
  std::cout << std::endl;
  std::cout << "Solution        = [";
  std::cout << finalPosition[0] << ",";
  std::cout << finalPosition[1] << "]" << std::endl;

  //
  // check results to see if it is within range
  //
  bool pass = true;
  double trueParameters[2] = { 1, 1 };
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      pass = false;
    }

  std::cout << "MaximumNumberOfIterations: " << itkOptimizer->GetNumberOfIterations();
  std::cout << std::endl;

  itkOptimizer->Print( std::cout );
  std::cout << "Stop description = "
            << itkOptimizer->GetStopConditionDescription() << std::endl << std::endl;

  if( !pass )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;


}
