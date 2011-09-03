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
#include "itkDisplacementFieldTransform.h"
#include "itkTranslationTransform.h"
#include "vnl/vnl_math.h"

/*
 * This test creates synthetic images and verifies numerical results
 * of metric evaluation.
 *
 * TODO
 * Test assigning displacement field of wrong size, expect exception.
 * Test with displacement field for fixed image transform.
 * Test evaluating over sub-region, maybe with non-identity tx's.
 * Test assigning different virtual image.
 * Test various options for image gradient calculation
 * Test permutations of image pre-warping and image gradient calculation method.
 * Test mask
 * Test with non-identity transforms
 * Test with gradient calculation performed in test itself rather than relying
 *  on the metric's gradient calculation.
 * Exercise other methods
 */

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
  typedef typename Superclass::FixedImageGradientType
                                                      FixedImageGradientType;
  typedef typename Superclass::MovingImagePointType   MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType   MovingImagePixelType;
  typedef typename Superclass::MovingImageGradientType
                                                     MovingImageGradientType;

  /* Implement pure virtual methods */
  void Initialize() throw ( itk::ExceptionObject )
  {
    //Be sure to call base class initialize
    Superclass::Initialize();

    //Now do your own initialization here
  }

  /* Provide the worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientType &  fixedImageGradient,
                    const MovingImagePointType &    mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientType & movingImageGradient,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType )
  {
    /* Just return some test values that can verify proper mechanics */
    metricValueResult = fixedImageValue + movingImageValue;
    
    //std::cout << "TestDerivedMetric: in GetValueAndDerivativeProcessPoint."
    //            << std::endl;
    if( mappedMovingPoint[0] > 47 && mappedMovingPoint[0] < 52 &&
        mappedMovingPoint[1] > 47 && mappedMovingPoint[1] < 52 )
        {
        std::cout << " mappedMovingPoint: " << mappedMovingPoint
                  << ": movingImageDerivative:"
                  << movingImageGradient << std::endl;
        }
        /*      << " mappedFixedPoint:  " << mappedFixedPoint
              << ":  fixedImageDerivative: "
              << fixedImageGradient << std::endl;
        */
    for ( unsigned int par = 0;
          par < this->GetNumberOfLocalParameters(); par++ )
      {
      double sum = 0.0;
      for ( unsigned int dim = 0;
            dim < Superclass::MovingImageDimension; dim++ )
        {
        sum += movingImageGradient[dim] + fixedImageGradient[dim];
        }
      localDerivativeReturn[par] = sum;
      }
    //  std::cout << " localDerivativeReturn: " << localDerivativeReturn << std::endl;

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
  void PrintSelf(std::ostream&, Indent) const {}

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
const unsigned int imageSize = 100;
const unsigned int imageDimensionality = 2;
typedef Image< double, imageDimensionality >              ImageType;
//typedef Transform< double, imageDimensionality, imageDimensionality >
//                                                          TransformType;
typedef TestDerivedMetric<ImageType,ImageType,ImageType>  TestMetricType;


//
// Compute truth values for the identity-transform tests
//
int ComputeIdentityTruthValues(const TestMetricType::Pointer & metric,
                               const ImageType::Pointer & fixedImage,
                               const ImageType::Pointer & movingImage,
                               TestMetricType::MeasureType& truthValue,
                               TestMetricType::DerivativeType& truthDerivative )
{
  // Make sure the metric is initialized
  metric->Initialize();
  // Call once to setup gradient images if applicable
  TestMetricType::MeasureType     tempValue;
  TestMetricType::DerivativeType  tempDerivative;
  
  metric->GetValueAndDerivative( tempValue, tempDerivative );

  truthValue = 0;
  truthDerivative.SetSize( metric->GetNumberOfParameters() );
  truthDerivative.Fill( 0 );

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

    // Get the image derivatives.
    // Because this test is using identity transforms,
    // simply retrieve by index.
    // NOTE: relying on the metric's gradient image isn't a
    // complete test, but it does test the rest of the mechanics.
    TestMetricType::MovingImageGradientType movingImageDerivative;
    TestMetricType::FixedImageGradientType  fixedImageDerivative;
    if( metric->GetUseFixedGradientRecursiveGaussianImageFilter() )
      {
      TestMetricType::FixedGradientImageType::ConstPointer
                   fixedGradientImage = metric->GetFixedGaussianGradientImage();
      fixedImageDerivative = fixedGradientImage->GetPixel( itFixed.GetIndex() );
      }
    else
      {
      typedef TestMetricType::FixedGradientCalculatorType::ConstPointer
                                              FixedGradientCalculatorPointer;
      FixedGradientCalculatorPointer fixedGradientCalculator =
                                        metric->GetFixedGradientCalculator();
      fixedImageDerivative =
              fixedGradientCalculator->EvaluateAtIndex( itFixed.GetIndex() );
      // We can skip the call to TransformCovariantVector since we're
      // working with identity transforms only.
      }
    if( metric->GetUseMovingGradientRecursiveGaussianImageFilter() )
      {
      TestMetricType::MovingGradientImageType::ConstPointer
              movingGradientImage = metric->GetMovingGaussianGradientImage();
      movingImageDerivative =
                        movingGradientImage->GetPixel( itMoving.GetIndex() );
      }
    else
      {
      typedef TestMetricType::MovingGradientCalculatorType::ConstPointer
                                              MovingGradientCalculatorPointer;
      MovingGradientCalculatorPointer     movingGradientCalculator;
      movingGradientCalculator = metric->GetMovingGradientCalculator();        
      movingImageDerivative =
              movingGradientCalculator->EvaluateAtIndex( itMoving.GetIndex() );
      }
      
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
    
} // computeTruthValues

////////////////////////////////////////////////////////////
//
// Pass true for 'setTruthValues' to have the results of the metric
// call set the return values fro truthValue and truthDerivative.
// Useful for establishing a relative truth for multiple runs.
// Otherwise, this will compare the results of calling the metric
// with truthValue and truthDerivative.
int RunSingleTest( const TestMetricType::Pointer & metric,
             const ImageType::Pointer & fixedImage,
             const ImageType::Pointer & movingImage,
             TestMetricType::MeasureType& truthValue,
             TestMetricType::DerivativeType& truthDerivative,
             int expecteNumberOfPoints,
             bool setTruthValues )
{
  std::cout << "Pre-warp image: fixed, moving: "
            << metric->GetPreWarpFixedImage() << ", "
            << metric->GetPreWarpMovingImage() << std::endl
            << "Use RecursiveGaussian filter for: fixed, moving: "
            << metric->GetUseFixedGradientRecursiveGaussianImageFilter()
            << ", "
            << metric->GetUseMovingGradientRecursiveGaussianImageFilter()
            << std::endl;                    

  // Initialize.
  //std::cout << "Initializing Metric." << std::endl;
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
  //std::cout << "Calling GetValueAndDerivative..." << std::endl;
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

  //Check number of threads and valid points
  std::cout << "Number of threads used: "
            << metric->GetNumberOfThreads() << std::endl;
  if( metric->GetNumberOfValidPoints() != ( expecteNumberOfPoints ) )
    {
    std::cout << "Expected number of valid points to be "
              << expecteNumberOfPoints
              << " but instead got " << metric->GetNumberOfValidPoints()
              << std::endl;
    return EXIT_FAILURE;
    }
  
  // Return or verify results
  int result = EXIT_SUCCESS;
  if( setTruthValues )
    {
    truthValue = valueReturn;
    truthDerivative = derivativeReturn;
    }
  else
    {
    // Verify results
    //std::cout << "Verifying results." << std::endl;
    if( vcl_fabs( truthValue - valueReturn ) > epsilon )
      {
      std::cout << "-FAILED- truthValue does not equal value: " << std::endl
                << "truthValue: " << truthValue << std::endl
                << "value: " << valueReturn << std::endl;
      result = EXIT_FAILURE;

      }
    if( ! testArray( truthDerivative, derivativeReturn ) )
      {
      std::cout << "-FAILED- truthDerivative does not equal derivatives:"
                << std::endl
                << "truthDerivative: " << truthDerivative << std::endl
                << "derivatives: " << derivativeReturn << std::endl;
      result = EXIT_FAILURE;
      }
    }
  return result;
}

}//namespace

////////////////////////////////////////////////////////////
int itkImageToImageObjectMetricTest(int, char ** const)
{
  bool origGlobalWarningValue = Object::GetGlobalWarningDisplay();
  Object::SetGlobalWarningDisplay( true );
  
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
    itFixed.Set( count * count );
    //itFixed.Set( 1.0 );
    count++;
    ++itFixed;
    }
  ImageRegionIteratorWithIndex<ImageType> itMoving( movingImage, region );
  itMoving.GoToBegin();
  count = 1;
  while( !itMoving.IsAtEnd() )
    {
    itMoving.Set( count * count / 2.0 );
    //itMoving.Set(1.0);
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
  // Tell the metric to compute image gradients for both fixed and moving.
  metric->SetGradientSource( TestMetricType::GRADIENT_SOURCE_BOTH );

  // Enable ITK debugging output
  metric->SetDebug( false );

  // Evaluate the metric, using identity transforms.
  // Test with different numbers of threads.
  // Run through all the permutations of pre-warping and image gradient
  // calculation method.
  TestMetricType::MeasureType     truthValue;
  TestMetricType::DerivativeType  truthDerivative;
  for( ThreadIdType numberOfThreads = 1; numberOfThreads < 5; numberOfThreads++ )
    {
    metric->SetNumberOfThreads( numberOfThreads );
    for( char useGaussianMoving = 1;
            useGaussianMoving >= 0; useGaussianMoving-- )
      {
        for( char useGaussianFixed = 1;
                useGaussianFixed >= 0; useGaussianFixed-- )
        {
        //Have to recompute new truth values for each permutation of
        // image gradient calculation options.
        bool computeNewTruthValues = true;
        for( char preWarpFixed = 1; preWarpFixed >= 0; preWarpFixed-- )
          {
          for( char preWarpMoving = 1; preWarpMoving >= 0; preWarpMoving-- )
            {
            metric->SetPreWarpFixedImage( preWarpFixed == 1 );
            metric->SetPreWarpMovingImage( preWarpMoving == 1);
            metric->SetUseFixedGradientRecursiveGaussianImageFilter(
                                                         useGaussianFixed == 1);
            metric->SetUseMovingGradientRecursiveGaussianImageFilter(
                                                         useGaussianMoving == 1 );
 /*           if( computeNewTruthValues )
              {
              ComputeIdentityTruthValues(metric, fixedImage, movingImage,
                                          truthValue, truthDerivative );
              }
            std::cout << "* Testing with identity transforms..."
                      << std::endl;                    
            if( RunSingleTest( metric, fixedImage, movingImage,
                                truthValue, truthDerivative,
                                imageSize * imageSize, false )
                                                              != EXIT_SUCCESS )
              {
              result = EXIT_FAILURE;
              }
 */           computeNewTruthValues = false;
            }
          }
        }
      } // loop through permutations
    }

  //
  // Test with non-identity transforms
  // This will compare for the same results between permutations of
  // image gradient calculation options, and not against a ground truth.
  
  FixedTransformType::ParametersType params( 2 );
  params[0] = 0;
  params[1] = 0;
  fixedTransform->SetParameters( params );
  params[0] = 1;
  params[1] = 0;
  movingTransform->SetParameters( params );
  
  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );
  for( ThreadIdType numberOfThreads = 1; numberOfThreads < 2; numberOfThreads++ )
    {
    metric->SetNumberOfThreads( numberOfThreads );
    for( char useGaussianMoving = 1;
            useGaussianMoving >= 1; useGaussianMoving-- )
      {
        for( char useGaussianFixed = 1;
                useGaussianFixed >= 1; useGaussianFixed-- )
        {
        //Have to recompute new truth values for each permutation of
        // image gradient calculation options.
        bool setNewTruthValues = true;
        for( char preWarpFixed = 1; preWarpFixed >= 1; preWarpFixed-- )
          {
          for( char preWarpMoving = 1; preWarpMoving >= 0; preWarpMoving-- )
            {
            metric->SetPreWarpFixedImage( preWarpFixed == 1 );
            metric->SetPreWarpMovingImage( preWarpMoving == 1);
            metric->SetUseFixedGradientRecursiveGaussianImageFilter(
                                                         useGaussianFixed == 1);
            metric->SetUseMovingGradientRecursiveGaussianImageFilter(
                                                         useGaussianMoving == 1 );
I am here - trouble with this test cuz it seems cuz of edge conditions changing the derivative vals between pre-warped and orig calcs. If I do a large image 100x100 and only look at center points, they seems idential between orig and pre-warped calcs, so that's probably the issue.
             std::cout << "* Testing with NON-IDENTITY transforms ..."
                      << std::endl;                    
            if( RunSingleTest( metric, fixedImage, movingImage,
                                truthValue, truthDerivative,
                                //( imageSize-2 ) * (imageSize - 2),
                                9900,
                                setNewTruthValues )
                                                              != EXIT_SUCCESS )
              {
              result = EXIT_FAILURE;
              }
            setNewTruthValues = false;
            }
          }
        }
      } // loop through permutations
    }
  
  
  //
  // Test with an identity displacement field transform for moving image
  //

  // Create a displacement field transform
  typedef itk::DisplacementFieldTransform<double, imageDimensionality>
                                                    DisplacementTransformType;
  DisplacementTransformType::Pointer displacementTransform =
      DisplacementTransformType::New();
  typedef DisplacementTransformType::DisplacementFieldType FieldType;
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
  DisplacementTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );
  // Assign to transform
  displacementTransform->SetDisplacementField( field );

  // Assign it to the metric
  metric->SetMovingTransform( displacementTransform );
  
  fixedTransform->SetIdentity();
  metric->SetFixedTransform( fixedTransform );
  
  metric->SetPreWarpFixedImage( true );
  metric->SetPreWarpMovingImage( true );
  metric->SetUseFixedGradientRecursiveGaussianImageFilter( true );
  metric->SetUseMovingGradientRecursiveGaussianImageFilter( true );
  // Tell the metric to compute image gradients for both fixed and moving.
  metric->SetGradientSource( TestMetricType::GRADIENT_SOURCE_BOTH );

  //Evaluate the metric
  std::cout
    << "* Testing with identity DisplacementFieldTransform for moving image..."
    << std::endl;
  ComputeIdentityTruthValues(metric, fixedImage, movingImage,
                               truthValue, truthDerivative );  
  if( RunSingleTest( metric, fixedImage, movingImage,
                      truthValue, truthDerivative,
                      imageSize * imageSize, false ) != EXIT_SUCCESS )
    {
    result = EXIT_FAILURE;
    }

  Object::SetGlobalWarningDisplay( origGlobalWarningValue );

  return result;
}
