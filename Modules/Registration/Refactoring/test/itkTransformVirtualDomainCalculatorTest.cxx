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

#include "itkFlipImageFilter.h"
#include "itkImageFileReader.h"
#include "itkTransformVirtualDomainCalculator.h"

template<unsigned int ImageDimension>
int TransformVirtualDomainCalculator( int itkNotUsed( argc ), char *argv[] )
{
  typedef itk::Image<char, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typedef itk::TransformVirtualDomainCalculator<ImageType, ImageType>
    DomainCalculatorType;
  typename DomainCalculatorType::Pointer domainCalculator =
    DomainCalculatorType::New();

  //
  // Test 1:  Get two input images and determine the virtual domain.
  //
  domainCalculator->SetInputImage1( reader1->GetOutput() );
  domainCalculator->SetInputImage2( reader2->GetOutput() );
  domainCalculator->SetUsePhysicalConsistency( false );

  try
    {
    domainCalculator->CalculateVirtualDomainParameters();
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << "Exception caught (test1) !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "------------------- Test1 --------------------" << std::endl;
  std::cout << "Origin 1: " << reader1->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Origin 2: " << reader2->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Direction 1: " << std::endl
    << reader1->GetOutput()->GetDirection() << std::flush;
  std::cout << "Direction 2: " << std::endl
    << reader2->GetOutput()->GetDirection() << std::flush;

  domainCalculator->Print( std::cout, 3 );

  //
  // Test 2:  Take the first image and set the direction to identity and origin
  // to all 0's.  Apply a known rotation to the second image and deterimine
  // if the virtual domain parameters are as expected without using
  // physical consistency.
  //
  typename DomainCalculatorType::OriginType origin;
  typename DomainCalculatorType::DirectionType direction;

  origin.Fill( 0 );
  reader1->GetOutput()->SetOrigin( origin );
  direction.SetIdentity();
  reader1->GetOutput()->SetDirection( direction );

  // Apply a rotation of 90 degrees around the z-axis;
  origin[0] = 5.0;
  reader2->GetOutput()->SetOrigin( origin );
  direction(0, 0) =  vcl_cos( 0.5 * vnl_math::pi );
  direction(0, 1) = -vcl_sin( 0.5 * vnl_math::pi );
  direction(1, 0) =  vcl_sin( 0.5 * vnl_math::pi );
  direction(1, 1) =  vcl_cos( 0.5 * vnl_math::pi );
  reader2->GetOutput()->SetDirection( direction );

  domainCalculator->SetInputImage1( reader1->GetOutput() );
  domainCalculator->SetInputImage2( reader2->GetOutput() );
  domainCalculator->SetUsePhysicalConsistency( false );

  try
    {
    domainCalculator->CalculateVirtualDomainParameters();
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << "Exception caught (test2) !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  origin[0] = 2.5;
  if( origin != domainCalculator->GetVirtualDomainOrigin() )
    {
    std::cerr << "Virtual domain origin is not what is expected." << std::endl;
    return EXIT_FAILURE;
    }
  direction(0, 0) =  0.5 * vcl_sqrt( 2.0 );
  direction(0, 1) = -0.5 * vcl_sqrt( 2.0 );
  direction(1, 0) =  0.5 * vcl_sqrt( 2.0 );
  direction(1, 1) =  0.5 * vcl_sqrt( 2.0 );
  typename DomainCalculatorType::DirectionType virtualDirection =
    domainCalculator->GetVirtualDomainDirection();
  if( ( ( direction - virtualDirection ).GetVnlMatrix().frobenius_norm() ) >
    1e-5 )
    {
    std::cerr << "Virtual domain direction is not what is expected." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "------------------- Test2 --------------------" << std::endl;
  std::cout << "Origin 1: " << reader1->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Origin 2: " << reader2->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Direction 1: " << std::endl
    << reader1->GetOutput()->GetDirection() << std::flush;
  std::cout << "Direction 2: " << std::endl
    << reader2->GetOutput()->GetDirection() << std::flush;

  domainCalculator->Print( std::cout, 3 );

  //
  // Test 3:  Take the first image and call FlipAxisImageFilter to change the
  // origin and direction matrix.  Turn on physical consistency and see if
  // the values are returned as expected.
  //

  origin.Fill( 0 );
  reader1->GetOutput()->SetOrigin( origin );
  direction.SetIdentity();
  reader1->GetOutput()->SetDirection( direction );

  origin[0] = 5.0;
  reader2->GetOutput()->SetOrigin( origin );
  reader2->GetOutput()->SetDirection( direction );

  typedef itk::FlipImageFilter<ImageType> FlipperType;
  typename FlipperType::Pointer flipper = FlipperType::New();
  typename FlipperType::FlipAxesArrayType flipAxes;
  flipAxes.Fill( 1 );
  flipper->SetInput( reader2->GetOutput() );
  flipper->SetFlipAxes( flipAxes );
  flipper->SetFlipAboutOrigin( false );
  flipper->Update();

  domainCalculator->SetInputImage1( reader1->GetOutput() );
  domainCalculator->SetInputImage2( flipper->GetOutput() );
  domainCalculator->SetUsePhysicalConsistency( true );

  // The virtual domain and origin should be exactly the same as for
  // Test2 since flip axes filter does not modify the actual physical location
  // of voxels in the image.

  try
    {
    domainCalculator->CalculateVirtualDomainParameters();
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << "Exception caught (test2) !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  origin[0] = 2.5;
  if( origin != domainCalculator->GetVirtualDomainOrigin() )
    {
    std::cerr << "Virtual domain origin is not what is expected." << std::endl;
    return EXIT_FAILURE;
    }
  direction.SetIdentity();
  virtualDirection = domainCalculator->GetVirtualDomainDirection();
  if( ( ( direction - virtualDirection ).GetVnlMatrix().frobenius_norm() ) >
    1e-5 )
    {
    std::cerr << "Virtual domain direction is not what is expected." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "------------------- Test3 --------------------" << std::endl;
  std::cout << "Origin 1: " << reader1->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Origin 2: " << flipper->GetOutput()->GetOrigin()
    << std::endl;
  std::cout << "Direction 1: " << std::endl
    << reader1->GetOutput()->GetDirection() << std::flush;
  std::cout << "Direction 2: " << std::endl
    << flipper->GetOutput()->GetDirection() << std::flush;

  domainCalculator->Print( std::cout, 3 );

  return EXIT_SUCCESS;
}

int itkTransformVirtualDomainCalculatorTest( int argc, char *argv[] )
{
  std::cout << "*******   " << argc << std::endl;

  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0]
      << " imageDimension inputImage1 inputImage2" << std::endl;
    exit( EXIT_FAILURE );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      return TransformVirtualDomainCalculator<2>( argc, argv );
      break;
    case 3:
      return TransformVirtualDomainCalculator<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
  return EXIT_SUCCESS;
}
