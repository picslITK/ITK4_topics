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
#include "itkDemonsImageToImageObjectMetric.h"
#include "itkGradientDescentObjectOptimizer.h"

#include "itkIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCommand.h"
#include "vnl/vnl_math.h"
#include "itkVectorCastImageFilter.h"

/* This test uses simple test images with circles created at runtime.
   Modified from itkDemonsRegistrationFilterTest */

namespace{
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
template<typename TRegistration>
class ShowProgressObject
{
public:
  ShowProgressObject(TRegistration* o)
    {m_Process = o;}
  void ShowProgress()
    {
    std::cout << "Progress: " << m_Process->GetProgress() << "  ";
    std::cout << "Iter: " << m_Process->GetElapsedIterations() << "  ";
    std::cout << "Metric: "   << m_Process->GetMetric()   << "  ";
    std::cout << "RMSChange: " << m_Process->GetRMSChange() << "  ";
    std::cout << std::endl;
    if ( m_Process->GetElapsedIterations() == 150 )
      { m_Process->StopRegistration(); }
    }
  typename TRegistration::Pointer m_Process;
};

// Template function to fill in an image with a circle.
template <class TImage>
void
FillWithCircle(
TImage * image,
double * center,
double radius,
typename TImage::PixelType foregnd,
typename TImage::PixelType backgnd )
{

  typedef itk::ImageRegionIteratorWithIndex<TImage> Iterator;
  Iterator it( image, image->GetBufferedRegion() );
  it.Begin();

  typename TImage::IndexType index;
  double r2 = vnl_math_sqr( radius );

  for( ; !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    double distance = 0;
    for( unsigned int j = 0; j < TImage::ImageDimension; j++ )
      {
      distance += vnl_math_sqr((double) index[j] - center[j]);
      }
    if( distance <= r2 ) it.Set( foregnd );
    else it.Set( backgnd );
    }

}

}//anonymous namespace

int itkDemonsRegistrationFilterTest(int, char* [] )
{

  typedef unsigned char PixelType;
  enum {ImageDimension = 2};
  typedef itk::Image<PixelType,ImageDimension> ImageType;
  typedef ImageType::IndexType  IndexType;
  typedef ImageType::SizeType   SizeType;
  typedef ImageType::RegionType RegionType;

  //--------------------------------------------------------
  std::cout << "Generate input images and deformation field";
  std::cout << std::endl;

  ImageType::SizeValueType sizeArray[ImageDimension] = { 128, 128 };
  SizeType size;
  size.SetSize( sizeArray );

  IndexType index;
  index.Fill( 0 );

  RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  ImageType::Pointer movingImage = ImageType::New();
  ImageType::Pointer fixedImage = ImageType::New();

  movingImage->SetLargestPossibleRegion( region );
  movingImage->SetBufferedRegion( region );
  movingImage->Allocate();

  fixedImage->SetLargestPossibleRegion( region );
  fixedImage->SetBufferedRegion( region );
  fixedImage->Allocate();

  double center[ImageDimension];
  double radius;
  PixelType fgnd = 250;
  PixelType bgnd = 15;

  // fill movingImage with circle
  center[0] = 64; center[1] = 64; radius = 30;
  FillWithCircle<ImageType>( movingImage, center, radius, fgnd, bgnd );

  // fill fixed with circle
  center[0] = 62; center[1] = 64; radius = 32;
  FillWithCircle<ImageType>( fixedImage, center, radius, fgnd, bgnd );

  //create a deformation field transform
  typedef DeformationFieldTransform<double, Dimension>
                                                    DeformationTransformType;
  DeformationTransformType::Pointer deformationTransform =
                                              DeformationTransformType::New();
  typedef DeformationTransformType::DeformationFieldType DeformationFieldType;
  DeformationFieldType::Pointer field = DeformationFieldType::New();

  //set the field to be the same as the fixed image region, which will
  // act by default as the virtual domain in this example.
  field->SetRegions( fixedImage->GetLargestPossibleRegion() );
  std::cout << "fixedImage->GetLargestPossibleRegion(): "
            << fixedImage->GetLargestPossibleRegion() << std::endl
            << "fixedImage->GetBufferedRegion(): "
            << fixedImage->GetBufferedRegion() << std::endl;
  field->Allocate();
  // Fill it with 0's
  DeformationTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );
  // Assign to transform
  deformationTransform->SetDeformationField( field );

  //identity transform for fixed image
  typedef IdentityTransform<double, Dimension> IdentityTransformType;
  IdentityTransformType::Pointer identityTransform =
                                                  IdentityTransformType::New();
  identityTransform->SetIdentity();

  /*
  typedef itk::VectorCastImageFilter<FieldType,FieldType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();

  //-------------------------------------------------------------
  std::cout << "Run registration and warp moving" << std::endl;

  typedef itk::DemonsRegistrationFilter<ImageType,ImageType,FieldType>
    RegistrationType;
  RegistrationType::Pointer registrator = RegistrationType::New();

  registrator->SetInitialDeformationField( caster->GetOutput() );
  registrator->SetMovingImage( movingImage );
  registrator->SetFixedImage( fixedImage );
  registrator->SetNumberOfIterations( 200 );
  registrator->SetStandardDeviations( 1.0 );
  registrator->SetMaximumError( 0.08 );
  registrator->SetMaximumKernelWidth( 10 );
  registrator->SetIntensityDifferenceThreshold( 0.001 );

  // turn on inplace execution
  registrator->InPlaceOn();

  // turn on/off use moving image gradient
  registrator->UseMovingImageGradientOff();

  typedef RegistrationType::DemonsRegistrationFunctionType FunctionType;
  FunctionType * fptr;
  fptr = dynamic_cast<FunctionType *>(
    registrator->GetDifferenceFunction().GetPointer() );
  fptr->Print( std::cout );

  // exercise other member variables
  std::cout << "No. Iterations: " << registrator->GetNumberOfIterations() << std::endl;
  std::cout << "Max. kernel error: " << registrator->GetMaximumError() << std::endl;
  std::cout << "Max. kernel width: " << registrator->GetMaximumKernelWidth() << std::endl;

  double v[ImageDimension];
  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    v[j] = registrator->GetStandardDeviations()[j];
    }
  registrator->SetStandardDeviations( v );
  */

am here in my adaptation...

  typedef ShowProgressObject<RegistrationType> ProgressType;
  ProgressType progressWatch(registrator);
  itk::SimpleMemberCommand<ProgressType>::Pointer command;
  command = itk::SimpleMemberCommand<ProgressType>::New();
  command->SetCallbackFunction(&progressWatch,
                               &ProgressType::ShowProgress);
  registrator->AddObserver( itk::ProgressEvent(), command);

  // warp moving image
  typedef itk::WarpImageFilter<ImageType,ImageType,FieldType> WarperType;
  WarperType::Pointer warper = WarperType::New();

  typedef WarperType::CoordRepType CoordRepType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,CoordRepType>
    InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();


  warper->SetInput( movingImage );
  warper->SetDeformationField( registrator->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixedImage->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOrigin() );
  warper->SetOutputDirection( fixedImage->GetDirection() );
  warper->SetEdgePaddingValue( bgnd );

  warper->Print( std::cout );

  warper->Update();

  // ---------------------------------------------------------
  std::cout << "Compare warped moving and fixed." << std::endl;

  // compare the warp and fixed images
  itk::ImageRegionIterator<ImageType> fixedIter( fixedImage,
      fixedImage->GetBufferedRegion() );
  itk::ImageRegionIterator<ImageType> warpedIter( warper->GetOutput(),
      fixedImage->GetBufferedRegion() );

  unsigned int numPixelsDifferent = 0;
  while( !fixedIter.IsAtEnd() )
    {
    if( fixedIter.Get() != warpedIter.Get() )
      {
      numPixelsDifferent++;
      }
    ++fixedIter;
    ++warpedIter;
    }

  std::cout << "Number of pixels different: " << numPixelsDifferent;
  std::cout << std::endl;

  if( numPixelsDifferent > 10 )
    {
    std::cout << "Test failed - too many pixels different." << std::endl;
    return EXIT_FAILURE;
    }

  registrator->Print( std::cout );

  // -----------------------------------------------------------
  std::cout << "Test running registrator without initial deformation field.";
  std::cout << std::endl;

  bool passed = true;
  try
    {
    registrator->SetInput( NULL );
    registrator->SetNumberOfIterations( 2 );
    registrator->Update();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    passed = false;
    }

  if ( !passed )
    {
    std::cout << "Test failed" << std::endl;
    return EXIT_FAILURE;
    }

  //--------------------------------------------------------------
  std::cout << "Test exception handling." << std::endl;

  std::cout << "Test NULL moving image. " << std::endl;
  passed = false;
  try
    {
    registrator->SetInput( caster->GetOutput() );
    registrator->SetMovingImage( NULL );
    registrator->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught expected error." << std::endl;
    std::cout << err << std::endl;
    passed = true;
    }

  if ( !passed )
    {
    std::cout << "Test failed" << std::endl;
    return EXIT_FAILURE;
    }
  registrator->SetMovingImage( movingImage );
  registrator->ResetPipeline();

  std::cout << "Test NULL moving image interpolator. " << std::endl;
  passed = false;
  try
    {
    fptr = dynamic_cast<FunctionType *>(
      registrator->GetDifferenceFunction().GetPointer() );
    fptr->SetMovingImageInterpolator( NULL );
    registrator->SetInput( initField );
    registrator->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught expected error." << std::endl;
    std::cout << err << std::endl;
    passed = true;
    }

  if ( !passed )
    {
    std::cout << "Test failed" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed" << std::endl;
  return EXIT_SUCCESS;


}
