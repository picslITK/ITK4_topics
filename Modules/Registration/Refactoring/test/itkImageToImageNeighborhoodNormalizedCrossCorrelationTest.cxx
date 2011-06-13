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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkVector.h"

// #include "itkObjectToObjectMetric.h"

#include "itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction.h"
#include "itkImageToData.h"


int main(int argc, char * argv[])
{

  if (argc <= 3) {
    std::cout << "Args: 2D_fixed_image 2D_moving_image radius_in_one_dim number_of_threads" << std::endl;
    exit(-1);
  }

  const char * filename1 = argv[1];
  const char * filename2 = argv[2];
  int radius_in_one_dim = atoi(argv[3]);
  int number_of_threads = atoi(argv[4]);





  const int ImageDimension = 2;
  typedef itk::Image< unsigned char, ImageDimension > ImageType;
  typedef ImageType::Pointer ImagePointerType;
  typedef ImageType::RegionType RegionType;

  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;


  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename1);
  reader->Update();
  ImagePointerType fixed_image = reader->GetOutput();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName(filename2);
  reader2->Update();
  ImagePointerType moving_image = reader2->GetOutput();

  VectorImageType::Pointer field = VectorImageType::New();
  field->SetRegions(fixed_image->GetLargestPossibleRegion());
  field->Allocate();

  VectorImageType::Pointer fieldInv = VectorImageType::New();
  fieldInv->SetRegions(fixed_image->GetLargestPossibleRegion());
  fieldInv->Allocate();

  typedef itk::ImageToImageNeighborhoodNormalizedCrossCorrelationFunction<ImageType, ImageType> ObjectMetricType;
  typedef itk::MetricThreadedHolder<ObjectMetricType, VectorImageType> MetricThreadedHolderType;
  typedef itk::ImageToData<ImageDimension, MetricThreadedHolderType> MetricThreaderType;


  itk::Size<ImageDimension> neighborhood_radius;
  neighborhood_radius.Fill(radius_in_one_dim);

  // pseudo code
  ObjectMetricType::Pointer objectMetric = ObjectMetricType::New();
  MetricThreadedHolderType metricHolder;
  MetricThreaderType::Pointer metricThreader = MetricThreaderType::New();


  objectMetric->SetFixedImage(fixed_image);
  objectMetric->SetMovingImage(moving_image);
  objectMetric->SetRadius(neighborhood_radius);
  objectMetric->InitializeGradientCalculator();

//  metricHolder.SetMetricFunction(objectMetric);
  metricHolder.metric = objectMetric;
  metricHolder.updateField = field;
  metricHolder.updateFieldInv = fieldInv;
  metricHolder.measure_per_thread.resize(number_of_threads);

  ImageType::RegionType inboundary_region = fixed_image->GetLargestPossibleRegion();

  metricThreader->SetNumberOfThreads(number_of_threads);
  metricThreader->m_OverallRegion = inboundary_region ;
  metricThreader->m_Holder = &metricHolder;
  metricThreader->ThreadedGenerateData = MetricThreadedHolderType::ComputeMetricValueInRegionOnTheFlyThreaded;

  metricThreader->GenerateData();

  float energy = static_cast<float> (metricHolder.AccumulateMeasuresFromAllThreads());

  std::cout << "cross correlation = " << energy << std::endl;

  return 1;

}
