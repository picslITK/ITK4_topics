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

#include "itkGradientDescentThreadedMetricOptimizer.h"
#include "itkDemonsImageToImageMetric.h"

using namespace itk;

/**
 */
int itkGradientDescentThreadedMetricOptimizerTest(int , char* [])
{
  const int ImageDimension = 2;
  typedef itk::Image<float, ImageDimension>                  ImageType;
  typedef ImageType::Pointer                                 ImagePointerType;
  typedef DemonsImageToImageMetric< ImageType, ImageType >   MetricType;
  typedef MetricType::RegionType                             ImageRegionType;

  MetricType::Pointer metric = MetricType::New();
  TestOptimizer::Pointer optimizer = TestOptimizer::New();
  ImageRegionType imageRegion;

  optimizer->SetMetric( metric );

  optimizer->SetNumberOfThreads(  );

  optimizer->SetOverallRegion( imageRegion );

  optimizer->StartOptimization();

  return EXIT_SUCCESS;
}
