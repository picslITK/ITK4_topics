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
#include "itkDemonsVectorImageToVectorImageObjectMetric.h"
#include "itkIdentityTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"

//These two are needed as long as we're using fwd-declarations in
//DisplacementFieldTransfor:
#include "itkVectorInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"

/* Simple test to verify that class builds and runs.
 * Results are not verified. See VectorImageToVectorImageObjectMetricTest
 * for verification of basic metric functionality.
 *
 * TODO Numerical verification.
 */

using namespace itk;

int itkDemonsVectorImageToVectorImageObjectMetricTest(int, char ** const)
{

  const unsigned int imageSize = 5;
  const unsigned int imageDimensionality = 3;
  typedef VectorImage< double, imageDimensionality >              ImageType;

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

  /* Create simple test images. */
  ImageType::Pointer fixedImage = ImageType::New();
  fixedImage->SetRegions( region );
  fixedImage->SetSpacing( spacing );
  fixedImage->SetOrigin( origin );
  fixedImage->SetDirection( direction );
  fixedImage->SetNumberOfComponentsPerPixel( 3 );
  fixedImage->Allocate();

  ImageType::Pointer movingImage = ImageType::New();
  movingImage->SetRegions( region );
  movingImage->SetSpacing( spacing );
  movingImage->SetOrigin( origin );
  movingImage->SetDirection( direction );
  movingImage->SetNumberOfComponentsPerPixel( 3 );
  movingImage->Allocate();

  /* Fill images */
  ImageRegionIterator<ImageType> itFixed( fixedImage, region );
  itFixed.GoToBegin();
  unsigned int count = 1;
  while( !itFixed.IsAtEnd() )
    {
    ImageType::PixelType fPix;
    fPix.SetSize(3);
    fPix.Fill( count*count );
    itFixed.Set( fPix );
    count++;
    ++itFixed;
    }
  ImageRegionIteratorWithIndex<ImageType> itMoving( movingImage, region );
  itMoving.GoToBegin();
  count = 1;
  while( !itMoving.IsAtEnd() )
    {
    ImageType::PixelType mPix;
    mPix.SetSize(3);
    mPix.Fill( 1.0/(count*count) );
    itMoving.Set( mPix );
    count++;
    ++itMoving;
    }

  /* Transforms */
  typedef TranslationTransform<double,imageDimensionality> FixedTransformType;
  typedef TranslationTransform<double,imageDimensionality> MovingTransformType;
  FixedTransformType::Pointer fixedTransform = FixedTransformType::New();
  MovingTransformType::Pointer movingTransform = MovingTransformType::New();
  fixedTransform->SetIdentity();
  movingTransform->SetIdentity();

  /* The metric */
  typedef DemonsVectorImageToVectorImageObjectMetric< ImageType, ImageType, ImageType >
                                                                  MetricType;
  MetricType::Pointer metric = MetricType::New();

  /* Assign images and transforms.
   * By not setting a virtual domain image or virtual domain settings,
   * the metric will use the fixed image for the virtual domain. */
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );

  /* Initialize. */
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
  //metric->SetNumberOfThreads( 1 );
  std::cout << "Threads: " << metric->GetNumberOfThreads() << std::endl;
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

  std::cout << "Value: " << valueReturn << std::endl;
  std::cout << "Derivative: " << derivativeReturn << std::endl;

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;
}
