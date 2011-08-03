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
#include "itkVectorCentralDifferenceImageFunction.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkVector.h"
#include "itkVectorImage.h"

int itkVectorCentralDifferenceImageFunctionTest(int, char* [] )
{
  const unsigned int ImageDimension = 2;
  typedef itk::Vector<unsigned int,3>          PixelType;
  typedef itk::Image<PixelType,ImageDimension> ImageType;

  typedef itk::VectorImage<unsigned int,2> VectorImageType;

  ImageType::Pointer image = ImageType::New();
  ImageType::SizeType size;
  size.Fill( 16 );
  ImageType::RegionType region( size );
  image->SetRegions( region );
  image->Allocate();

  VectorImageType::Pointer vimage = VectorImageType::New();
  VectorImageType::SizeType vsize;
  vsize.Fill( 16 );
  VectorImageType::RegionType vregion( vsize );
  vimage->SetRegions( vregion );
  vimage->SetNumberOfComponentsPerPixel( 3 );
  vimage->Allocate();

  // make a test image
  typedef itk::ImageRegionIterator<ImageType> Iterator;
  Iterator iter( image, region );
  iter.GoToBegin();
  unsigned int counter = 0;

  while ( !iter.IsAtEnd() )
    {
    PixelType pix;
    pix.Fill( counter );
    iter.Set( pix );
    ++counter;
    ++iter;
    }

  typedef itk::ImageRegionIterator<VectorImageType> VIterator;
  VIterator viter( vimage, vregion );
  viter.GoToBegin();
  counter = 0;

  while ( !viter.IsAtEnd() )
    {
    VectorImageType::PixelType vpix;
    vpix.SetSize( 3 );
    vpix.Fill( counter );
    viter.Set( vpix );
    ++counter;
    ++viter;
    }

  // set up central difference calculator
  typedef float CoordRepType;
  typedef itk::VectorCentralDifferenceImageFunction<ImageType,CoordRepType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();

  function->SetInputImage( image );

  ImageType::IndexType index;

  // pick an index inside the image
  index.Fill( 8 );
  std::cout << "Index: " << index << " Derivative: ";
  std::cout << function->EvaluateAtIndex( index ) << std::endl;

  if ( function->IsInsideBuffer( index ) )
    {
    std::cout << "Index: " << index << " is inside the BufferedRegion." << std::endl;
    }

  FunctionType::ContinuousIndexType cindex;
  cindex.Fill( 8.0 );
  std::cout << "ContinuousIndex: " << cindex << " Derivative: ";
  std::cout << function->EvaluateAtContinuousIndex( cindex ) << std::endl;

  FunctionType::PointType point;
  point.Fill( 8.0 );
  std::cout << "Point: " << cindex << " Derivative: ";
  std::cout << function->Evaluate( point ) << std::endl;

  // pick an index on the image edge
  index.Fill( 8 );
  index[0] = 15;
  std::cout << "Index: " << index << " Derivative: ";
  std::cout << function->EvaluateAtIndex( index ) << std::endl;

  if ( function->IsInsideBuffer( index ) )
    {
    std::cout << "Index: " << index << " is inside the BufferedRegion." << std::endl;
    }

  cindex.Fill( 8.0 );
  cindex[0] = 15.0;
  std::cout << "ContinuousIndex: " << cindex << " Derivative: ";
  std::cout << function->EvaluateAtContinuousIndex( cindex ) << std::endl;

  point.Fill( 8.0 );
  point[0] = 15.0;
  std::cout << "Point: " << cindex << " Derivative: ";
  std::cout << function->Evaluate( point ) << std::endl;


  // set up central difference calculator
  typedef float CoordRepType;
  typedef itk::VectorCentralDifferenceImageFunction<VectorImageType,CoordRepType> VectorFunctionType;
  VectorFunctionType::Pointer vfunction = VectorFunctionType::New();

  vfunction->SetInputImage( vimage );

  VectorImageType::IndexType vindex;

  // pick an index inside the image
  vindex.Fill( 8 );
  std::cout << "Index: " << vindex << " Derivative: ";
  std::cout << vfunction->EvaluateAtIndex( vindex ) << std::endl;

  if ( vfunction->IsInsideBuffer( vindex ) )
    {
    std::cout << "Index: " << vindex << " is inside the BufferedRegion." << std::endl;
    }

  VectorFunctionType::ContinuousIndexType vcindex;
  vcindex.Fill( 8.0 );
  std::cout << "ContinuousIndex: " << vcindex << " Derivative: ";
  std::cout << vfunction->EvaluateAtContinuousIndex( vcindex ) << std::endl;

  VectorFunctionType::PointType vpoint;
  vpoint.Fill( 8.0 );
  std::cout << "Point: " << vcindex << " Derivative: ";
  std::cout << vfunction->Evaluate( vpoint ) << std::endl;

  // pick an index on the image edge
  vindex.Fill( 8 );
  vindex[0] = 15;
  std::cout << "Index: " << vindex << " Derivative: ";
  std::cout << vfunction->EvaluateAtIndex( vindex ) << std::endl;

  if ( vfunction->IsInsideBuffer( vindex ) )
    {
    std::cout << "Index: " << vindex << " is inside the BufferedRegion." << std::endl;
    }

  vcindex.Fill( 8.0 );
  vcindex[0] = 15.0;
  std::cout << "ContinuousIndex: " << vcindex << " Derivative: ";
  std::cout << vfunction->EvaluateAtContinuousIndex( vcindex ) << std::endl;

  vpoint.Fill( 8.0 );
  vpoint[0] = 15.0;
  std::cout << "Point: " << vcindex << " Derivative: ";
  std::cout << vfunction->Evaluate( vpoint ) << std::endl;

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;

}
