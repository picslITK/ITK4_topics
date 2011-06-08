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

#include "itkObjectToObjectThreadedMetricOptimizerBase.h"
#include "itkDemonsImageToImageMetric.h"

using namespace itk;

namespace{

const int ImageDimension = 2;
typedef itk::Image<float, ImageDimension>                     ImageType;
typedef ImageType::Pointer                                    ImagePointerType;
typedef DemonsImageToImageMetric< ImageType, ImageType >      MetricType;
typedef MetricType::RegionType                                ImageRegionType;

/* Define a simple derived class. */
class TestOptimizer
  : public ObjectToObjectThreadedMetricOptimizerBase< MetricType >
{
public:
  /** Standard "Self" typedef. */
  typedef TestOptimizer                           Self;
  typedef ObjectToObjectThreadedMetricOptimizerBase< MetricType >
                                                  Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /* Provide an override for the pure virtual StartOptimization */
  void StartOptimization()
    {
    std::cout << "StartOptimization called." << std::endl;
    }

};

}//namespace

/**
 */
int itkObjectToObjectThreadedMetricOptimizerBaseTest(int , char* [])
{
  MetricType::Pointer metric = MetricType::New();
  TestOptimizer::Pointer optimizer = TestOptimizer::New();
  ImageRegionType imageRegion;

  /* exercise some methods */
  optimizer->SetMetric( metric );
  MetricType::Pointer metricReturn = optimizer->GetMetric();
  if( metricReturn != metric )
    {
    std::cerr << "Set/GetMetric failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "value: " << optimizer->GetValue() << std::endl;

  TestOptimizer::MetricThreaderPointer threader =
    optimizer->GetMetricThreader();
  /*check that the threader's holder has been properly set */
  if( threader->GetHolder() != static_cast<void*>(optimizer.GetPointer()) )
    {
    std::cerr << "Set/GetHolder failed." << std::endl;
    std::cerr << "GetHolder(): " << threader->GetHolder() << std::endl
              << "optimizer.GetPointer(): " << optimizer.GetPointer()
              << std::endl;
    return EXIT_FAILURE;
    }

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

  optimizer->SetOverallRegion( imageRegion );

  optimizer->StartOptimization();

  return EXIT_SUCCESS;
}
