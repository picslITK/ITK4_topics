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

#include "itkObjectToObjectOptimizerBase.h"
#include "itkDemonsImageToImageObjectMetric.h"

using namespace itk;

namespace{

const int ImageDimension = 2;
typedef itk::Image<double, ImageDimension>                    ImageType;
typedef ImageType::Pointer                                    ImagePointerType;

/* Define a simple derived class. */
class TestOptimizer
  : public ObjectToObjectOptimizerBase
{
public:
  /** Standard "Self" typedef. */
  typedef TestOptimizer                           Self;
  typedef ObjectToObjectOptimizerBase             Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TestOptimizer, ObjectToObjectOptimizerBase);

  /* Provide an override for the pure virtual StartOptimization */
  void StartOptimization()
    {
    std::cout << "StartOptimization called." << std::endl;
    }

};

}//namespace

/**
 */
int itkObjectToObjectOptimizerBaseTest(int , char* [])
{
  typedef DemonsImageToImageObjectMetric<ImageType,ImageType>   MetricType;
  MetricType::Pointer metric = MetricType::New();
  TestOptimizer::Pointer optimizer = TestOptimizer::New();

  /* exercise some methods */
  optimizer->SetMetric( metric );
  if( optimizer->GetMetric() != metric )
    {
    std::cerr << "Set/GetMetric failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "value: " << optimizer->GetValue() << std::endl;

  /* Test set/get of scales */
  TestOptimizer::ScalesType scales;
  scales.Fill(3);
  optimizer->SetScales( scales );
  const TestOptimizer::ScalesType& scalesReturn = optimizer->GetScales();
  if( scalesReturn != scales )
    {
    std::cerr << "Set/GetScales failed." << std::endl;
    return EXIT_FAILURE;
    }

  optimizer->SetNumberOfThreads( 1 );

  optimizer->StartOptimization();

  return EXIT_SUCCESS;
}
