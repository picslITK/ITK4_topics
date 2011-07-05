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
#include "itkImage.h"
#include "itkImageToImageObjectMetric.h"
#include "itkDeformationFieldTransform.h"
#include "itkTranslationTransform.h"
#include "vnl/vnl_math.h"

using namespace itk;

namespace {

template<class TFixedImage,class TMovingImage,class TVirtualImage>
class TestDerivedMetric
  : public ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:
  /** Standard class typedefs. */
  typedef TestDerivedMetric                                   Self;
  typedef ImageToImageObjectMetric<TFixedImage, TMovingImage,
                                               TVirtualImage> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TestDerivedMetric, ImageToImageObjectMetric);

  /** superclass types */
  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::VirtualPointType           VirtualPointType;
  typedef typename Superclass::FixedImagePointType        FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType        FixedImagePixelType;
  typedef typename Superclass::FixedImageDerivativesType
                                                     FixedImageDerivativesType;
  typedef typename Superclass::MovingImagePointType       MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType       MovingImagePixelType;
  typedef typename Superclass::MovingImageDerivativesType
                                                    MovingImageDerivativesType;

  /* Implement pure virtual methods */
  void Initialize() throw ( itk::ExceptionObject )
  {
    //Be sure to call base class initialize
    Superclass::Initialize();

    //Now do your own initialization here
  }

  /* Provide the worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  fixedImageDerivatives,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID)
  {
    /* Just return some test values that can verify proper mechanics */
    metricValueResult = fixedImageValue + movingImageValue;

    /*
    std::cout << "Metric: " << mappedMovingPoint << ": movingImageDerivative:"
              << movingImageDerivatives << std::endl
              << "Metric: " << mappedFixedPoint << ": fixedImageDerivative: "
              << fixedImageDerivatives << std::endl;
    */

    for ( unsigned int par = 0;
          par < this->GetNumberOfLocalParameters(); par++ )
      {
      double sum = 0.0;
      for ( unsigned int dim = 0;
            dim < Superclass::MovingImageDimension; dim++ )
        {
        sum += movingImageDerivatives[dim] + fixedImageDerivatives[dim];
        }
      localDerivativeReturn[par] = sum;
      }
      /*
      std::cout << "TestDerivedMetric: in GetValueAndDerivativeProcessPoint."
                << std::endl
                << "localDerivativeReturn: " << localDerivativeReturn << std::endl;
      */
    // Return true if the point was used in evaluation
    return true;
  }

  //This is of one two evaluation methods that the user may call.
  MeasureType GetValue()
  {
    //TODO
    return 0.0;
  }

  //This is of one two evaluation methods that the user may call.
  void GetValueAndDerivative( MeasureType & valueReturn,
                              DerivativeType & derivativeReturn)
  {
    //1) Do any pre-processing required for your metric. To help with
    // threading, you can use ImageToData or Array1DToData classes,
    // or derive your own from ObjectToData.

    //2) Call GetValueAndDerivativeMultiThreadedInitiate.
    //This will iterate over virtual image region and call your
    // GetValueAndDerivativeProcessPoint method, see definition in
    // base.
    this->GetValueAndDerivativeMultiThreadedInitiate( derivativeReturn );

    //3) Optionally call GetValueAndDerivativeMultiThreadedPostProcess for
    // default post-processing, which sums up results from each thread,
    // and optionally averages them. It then assigns the results to
    // 'value' and 'derivative', without copying in the case of 'derivative'.
    //Do your own post-processing as needed.
    this->GetValueAndDerivativeMultiThreadedPostProcess( true /*doAverage*/ );

    //4) Return the value result. The derivative result has already been
    // written to derivativeReturn.
    valueReturn = this->GetValueResult();

    //That's it. Easy as 1, 2, 3 (and 4).
  }

protected:
  TestDerivedMetric(){};
  virtual ~TestDerivedMetric() {}
  void PrintSelf(std::ostream& os, Indent indent) const {}

private:
  //purposely not implemented
  TestDerivedMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

}; //TestDerived Metric ///////////////////////////////////////////////////

const double epsilon = 1e-10;
template <typename TVector>
bool testArray( const TVector & v1, const TVector & v2 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < v1.Size(); i++ )
    {
    if( vcl_fabs( v1[i] - v2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }


//Global types
const unsigned int imageSize = 4;
const unsigned int imageDimensionality = 2;
typedef Image< double, imageDimensionality >              ImageType;
//typedef Transform< double, imageDimensionality, imageDimensionality >
//                                                          TransformType;
typedef TestDerivedMetric<ImageType,ImageType,ImageType>  TestMetricType;

////////////////////////////////////////////////////////////
int RunTest( TestMetricType::Pointer & metric,
             ImageType::Pointer & fixedImage,
             ImageType::Pointer & movingImage )
{
  // Initialize.
  try
    {
    metric->Initialize();
    }
  catch( ExceptionObject & exc )
    {
    std::cout << "Caught unexpected exception during Initialize: "
              << exc;
    return EXIT_FAILURE;
    }

  // Evaluate
  TestMetricType::MeasureType valueReturn;
  TestMetricType::DerivativeType derivativeReturn;
  std::cout << "Calling GetValueAndDerivative..." << std::endl;
  try
    {
    metric->GetValueAndDerivative( valueReturn, derivativeReturn );
    }
  catch( ExceptionObject & exc )
    {
    std::cout << "Caught unexpected exception during GetValueAndDerivative: "
              << exc;
    return EXIT_FAILURE;
    }
  std::cout << "...done GetValueAndDerivative." << std::endl;

  //Check number of threads and valid points
  std::cout << "Number of threads used: "
            << metric->GetNumberOfThreads() << std::endl;
  if( metric->GetNumberOfValidPoints() != (imageSize * imageSize ) )
    {
    std::cout << "Expected number of valid points to be "
              << imageSize * imageSize
              << " but instead got " << metric->GetNumberOfValidPoints();
    return EXIT_FAILURE;
    }

  //
  // Verify results
  //
  TestMetricType::MeasureType       truthValue = 0;
  TestMetricType::DerivativeType
                    truthDerivative( metric->GetNumberOfParameters() );
  truthDerivative.Fill( 0 );

  // Default behavior is for the metric to precompute image derivatives, so
  // we access them here for testing.
  TestMetricType::FixedImageDerivativesType  fixedImageDerivative;
  TestMetricType::MovingImageDerivativesType movingImageDerivative;
  TestMetricType::FixedGradientImageType::ConstPointer
    fixedGradientImage = metric->GetFixedGaussianGradientImage();
  TestMetricType::MovingGradientImageType::ConstPointer
    movingGradientImage = metric->GetMovingGaussianGradientImage();

  // Compute truth values
  ImageRegionIterator<ImageType>
    itFixed( fixedImage, fixedImage->GetRequestedRegion() );
  ImageRegionIterator<ImageType>
    itMoving( movingImage, movingImage->GetRequestedRegion() );
  itFixed.GoToBegin();
  itMoving.GoToBegin();
  unsigned int count = 0;
  while( !itFixed.IsAtEnd() && !itMoving.IsAtEnd() )
    {
    truthValue += itFixed.Get() + itMoving.Get();

    // Get the image derivatives. Because this test is using identity transforms,
    // simply retrieve by index.
    fixedImageDerivative = fixedGradientImage->GetPixel( itFixed.GetIndex() );
    movingImageDerivative = movingGradientImage->GetPixel( itMoving.GetIndex() );
    /*
    std::cout << "Truth: " << itMoving.GetIndex() << ": movingImageDerivative:"
              << movingImageDerivative << std::endl
              << "Truth: " << itFixed.GetIndex() << ": fixedImageDerivative: "
              << fixedImageDerivative << std::endl;
    */
    for ( unsigned int par = 0;
          par < metric->GetNumberOfLocalParameters(); par++ )
      {
      double sum = 0.0;
      for ( unsigned int dim = 0; dim < imageDimensionality; dim++ )
        {
        sum += movingImageDerivative[dim] + fixedImageDerivative[dim];
        }

      if( metric->HasLocalSupport() )
        {
        truthDerivative[ count * metric->GetNumberOfLocalParameters() + par ]
                                                                        = sum;
        }
      else
        {
        truthDerivative[par] += sum;
        }
      }
    count++;
    ++itFixed;
    ++itMoving;
    }

  // Take the averages
  truthValue /= metric->GetNumberOfValidPoints();
  if( ! metric->HasLocalSupport() )
    {
    truthDerivative /= metric->GetNumberOfValidPoints();
    }

  // Compare truth with test
  if( vcl_fabs( truthValue - valueReturn ) > epsilon )
    {
    std::cout << "truthValue does not equal value: " << std::endl
              << "truthValue: " << truthValue << std::endl
              << "value: " << valueReturn << std::endl;
    return EXIT_FAILURE;

    }
  if( ! testArray( truthDerivative, derivativeReturn ) )
    {
    std::cout << "truthDerivative does not equal derivatives:" << std::endl
              << "truthDerivative: " << truthDerivative << std::endl
              << "derivatives: " << derivativeReturn << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

}//namespace

////////////////////////////////////////////////////////////
int itkImageToImageObjectMetricTest(int argc, char * argv[])
{
  int result = EXIT_SUCCESS;

  ImageType::SizeType       size = {{imageSize, imageSize}};
  ImageType::IndexType      index = {{0,0}};
  ImageType::RegionType     region;
  region.SetSize( size );
  region.SetIndex( index );
  ImageType::SpacingType    spacing;
  spacing.Fill(1.0);
  ImageType::PointType      origin;
  origin.Fill(0);
  ImageType::DirectionType  direction;
  direction.SetIdentity();

  // Create simple test images.
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
  ImageRegionIterator<ImageType> itFixed( fixedImage, region );
  itFixed.GoToBegin();
  unsigned int count = 1;
  while( !itFixed.IsAtEnd() )
    {
    itFixed.Set( count*count );
    count++;
    ++itFixed;
    }
  ImageRegionIteratorWithIndex<ImageType> itMoving( movingImage, region );
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
  FixedTransformType::Pointer fixedTransform = FixedTransformType::New();
  MovingTransformType::Pointer movingTransform = MovingTransformType::New();
  fixedTransform->SetIdentity();
  movingTransform->SetIdentity();

  // The simplistic test metric
  TestMetricType::Pointer metric = TestMetricType::New();

  // Assign images and transforms.
  // By not setting a virtual domain image or virtual domain settings,
  // the metric will use the fixed image for the virtual domain.
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );

  //Evaluate the metric
  //metric->SetNumberOfThreads(1);
  std::cout << "==Testing with IdentityTransform for moving image..." << std::endl;
  if( RunTest( metric, fixedImage, movingImage ) != EXIT_SUCCESS )
    {
    result = EXIT_FAILURE;
    }

  //
  // Test with an identity deformation field transform for moving image
  //

  // Create a deformation field transform
  typedef itk::DeformationFieldTransform<double, imageDimensionality>
                                                    DeformationTransformType;
  DeformationTransformType::Pointer deformationTransform =
      DeformationTransformType::New();
  typedef DeformationTransformType::DeformationFieldType FieldType;
  FieldType::Pointer field = FieldType::New(); //This is based on itk::Image

  FieldType::SizeType defsize;
  FieldType::IndexType start;
  FieldType::RegionType defregion;
  defsize.Fill( imageSize );
  start.Fill( 0 );
  defregion.SetSize( defsize );
  defregion.SetIndex( start );
  field->SetRegions( defregion );
  field->Allocate();
  // Fill it with 0's
  DeformationTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );
  // Assign to transform
  deformationTransform->SetDeformationField( field );

  // Assign it to the metric
  metric->SetMovingTransform( deformationTransform );
  //Evaluate the metric
  std::cout
    << "==Testing with identity DeformationFieldTransform for moving image..."
    << std::endl;
  if( RunTest( metric, fixedImage, movingImage ) != EXIT_SUCCESS )
    {
    result = EXIT_FAILURE;
    }


  /////////////// TODO

  // Test assigning deformation field of wrong size, expect exception

  // Test with deformation field for fixed image tx

  // Test evaluating over sub-region, maybe with non-identity tx's

  // Test assigning different virtual image

  // Exercise other methods

  return result;
}
