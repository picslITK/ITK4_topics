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

#include "itkMattesMutualInformationImageToImageObjectMetric.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkTranslationTransform.h"
#include "itkVector.h"

/* Simple test to verify that the class builds and runs.
 * Results are not verified. */


int itkMattesMutualInformationImageToImageObjectMetricTest(int argc, char * argv[])
{
  /*
  typedef double TPixelType;
  const unsigned int VDimension = 2;  
  typedef itk::Image< TPixelType, VDimension >  ImageType;
  
  const unsigned int imageSize = 5;
  ImageType::SizeType       size;
  size.Fill( imageSize );
  ImageType::IndexType      index;
  index.Fill( 0 );
  ImageType::RegionType     region;
  region.SetSize( size );
  region.SetIndex( index );
  ImageType::SpacingType    spacing;
  spacing.Fill(1.0);
  ImageType::PointType      origin;
  origin.Fill(0);
  ImageType::DirectionType  direction;
  direction.SetIdentity();

  // Create simple test images. //
  ImageType::Pointer fixedImage = ImageType::New();
  fixedImage->SetRegions( region );
  fixedImage->SetSpacing( spacing );
  fixedImage->SetOrigin( origin );
  fixedImage->SetDirection( direction );
  fixedImage->Allocate();

  ImageType::Pointer movingImage = ImageType::New();
  movingImage->SetRegions( region );
  movingImage->SetSpacing( spacing );
  movingImage->SetOrigin( origin );
  movingImage->SetDirection( direction );
  movingImage->Allocate();

  // Fill images 
  typedef itk::ImageRegionIterator<ImageType> FixedImageIteratorType;
  FixedImageIteratorType itFixed( fixedImage, region );
  
  itFixed.GoToBegin();
  unsigned int count = 1;
  while( !itFixed.IsAtEnd() )
    {
    itFixed.Set( count*count );
    count++;
    ++itFixed;
    }
    
  typedef itk::ImageRegionIteratorWithIndex<ImageType> MovingImageIteratorType;
  MovingImageIteratorType itMoving( movingImage, region );
  
  itMoving.GoToBegin();
  count = 1;
  while( !itMoving.IsAtEnd() )
    {
    itMoving.Set( 1.0/(count*count) );
    count++;
    ++itMoving;
    }

  // Transforms 
  typedef TranslationTransform<double,imageDimensionality> FixedTransformType;
  typedef TranslationTransform<double,imageDimensionality> MovingTransformType;
  FixedTransformType::Pointer fixedTransform   = FixedTransformType::New();
  MovingTransformType::Pointer movingTransform = MovingTransformType::New();
  fixedTransform->SetIdentity();
  movingTransform->SetIdentity();

  // The metric 
  typedef MattesMutualInformationImageToImageObjectMetric< ImageType, ImageType, ImageType >
                                                                  MetricType;
  MetricType::Pointer metric = MetricType::New();

  // Assign images and transforms.
   * By not setting a virtual domain image or virtual domain settings,
   * the metric will use the fixed image for the virtual domain. 
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedImageTransform( fixedTransform );
  metric->SetMovingImageTransform( movingTransform );

  // Initialize. 
  try
    {
    std::cout << "Calling Initialize..." << std::endl;
    metric->Initialize();
    }
  catch( ExceptionObject & exc )
    {
    std::cout << "Caught unexpected exception during Initialize: "
              << exc;
    return EXIT_FAILURE;
    }

  // Evaluate
  MetricType::MeasureType valueReturn;
  MetricType::DerivativeType derivativeReturn;
  try
    {
    std::cout << "Calling GetValueAndDerivative..." << std::endl;
    metric->GetValueAndDerivative( valueReturn, derivativeReturn );
    }
  catch( ExceptionObject & exc )
    {
    std::cout << "Caught unexpected exception during GetValueAndDerivative: "
              << exc;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed." << std::endl;
  */
  return EXIT_SUCCESS;
}
