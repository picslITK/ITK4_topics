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

#include <iostream>
#include "itkVariableLengthVector.h"
#include "itkFixedArray.h"
#include "itkTimeProbe.h"
#include "itkVectorImageToImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"


// This test tests:
// VectorImage
// -----------
//  - Iterators on the VectorImage,
//  - Timing tests comparing itk::VectorImage to similar itk::Image using
//    FixedArray and VariableLengthVector
//  - IO support for VectorImage.

int itkVectorImageTest( int, char* argv[] )
{
  bool failed = false;

  const unsigned int Dimension    = 3;
  const unsigned int VectorLength = 2 * Dimension;
  typedef float PixelType;

  {
  // Test 1.
  //
  // Create an Image of VariableLengthVector, FixedArray, VectorImage of length 6 and compare
  // times.
  //
  // Three images.. for crude timing analysis.

  typedef itk::Image< itk::VariableLengthVector< PixelType >, Dimension > VariableLengthVectorImageType;
  typedef itk::Image< itk::FixedArray< PixelType, VectorLength >,
                                      Dimension > FixedArrayImageType;
  typedef itk::VectorImage< PixelType, Dimension >   VectorImageType;



  // Using image of VariableLengthVector< PixelType >
  {
  typedef itk::VariableLengthVector< PixelType > InternalPixelType;

  itk::TimeProbe clock;
  clock.Start();

  VariableLengthVectorImageType::Pointer image = VariableLengthVectorImageType::New();
  VariableLengthVectorImageType::IndexType start;
  InternalPixelType f( VectorLength );
  VariableLengthVectorImageType::SizeType  size;
  for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i; }
  start.Fill(0);
  size.Fill(50);
  VariableLengthVectorImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer( f );

  clock.Stop();
  double timeTaken = clock.GetMeanTime();
  std::cout << "Allocating an image of itk::VariableLengthVector of length " <<  VectorLength
          << " with image size " << size << " took " << timeTaken << " s." << std::endl;

    // Const iterator over the image...
    {
    clock.Start();
    typedef itk::ImageRegionConstIterator< VariableLengthVectorImageType > IteratorType;
    IteratorType it( image, image->GetBufferedRegion() );
    it.Begin();
    while( !it.IsAtEnd() )
      {
      it.Get();
      ++it;
      }
    clock.Stop();
    std::cout << "ConstIterator Get() over the entire image took : " <<
      clock.GetMeanTime() << " s." << std::endl;
    }
  }


  // Using image of FixedArray< PixelType, VectorLength >
  {
  itk::TimeProbe clock;
  clock.Start();

  typedef itk::FixedArray< PixelType, VectorLength > InternalPixelType;

  FixedArrayImageType::Pointer image = FixedArrayImageType::New();
  FixedArrayImageType::IndexType start;
  InternalPixelType f;
  FixedArrayImageType::SizeType  size;
  for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i; }
  start.Fill(0);
  size.Fill(50);
  FixedArrayImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer( f );

  clock.Stop();
  double timeTaken = clock.GetMeanTime();
  std::cout << "Allocating an image of itk::FixedArray of length " <<  VectorLength
          << " with image size " << size << " took " << timeTaken << " s." << std::endl;

    {
    // Test and compare times with iterators
    //
    // First set some pixel
    for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i*0.1; }
    FixedArrayImageType::IndexType idx;
    for (unsigned int i = 0; i < Dimension; i++)
      {
      idx[i] = 4;
      }
    idx[Dimension-1] = 12;
    image->SetPixel( idx, f );
    }

    {
    clock.Start();
    typedef itk::ImageRegionConstIterator< FixedArrayImageType > IteratorType;
    IteratorType it( image, image->GetBufferedRegion() );
    it.Begin();
    while( !it.IsAtEnd() )
      {
      it.Get();
      ++it;
      }
    clock.Stop();
    std::cout << "ConstIterator Get() over the entire image took : " <<
      clock.GetMeanTime() << " s." << std::endl;
    }
  }

  // Using VectorImage< PixelType, Dimension >
  {
  itk::TimeProbe clock;
  clock.Start();

  VectorImageType::Pointer vectorImage = VectorImageType::New();
  VectorImageType::IndexType start;
  itk::VariableLengthVector< PixelType > f( VectorLength );
  VectorImageType::SizeType  size;
  for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i; }
  start.Fill(0);
  size.Fill(50);
  VectorImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  vectorImage->SetVectorLength( VectorLength );
  vectorImage->SetRegions( region );
  vectorImage->Allocate();
  vectorImage->FillBuffer( f );

  clock.Stop();
  double timeTaken = clock.GetMeanTime();
  std::cout << "Allocating an image of itk::VectorImage with pixels length " <<  VectorLength
     << " with image size " << size << " took " << timeTaken << " s." << std::endl;


  // Iterator tests on the vector image.
  //
  // Const iterator over the vector image...
  {
    clock.Start();
    typedef itk::ImageRegionConstIterator< VectorImageType > IteratorType;
    IteratorType it( vectorImage, vectorImage->GetBufferedRegion() );
    it.Begin();
    while( !it.IsAtEnd() )
      {
      it.Get();
      ++it;
      }
    clock.Stop();
    std::cout << "ConstIterator Get() over the entire vectorImage took : " <<
      clock.GetMeanTime() << " s." << std::endl;
  }



  std::cout << "---------------------------------------------------------------" << std::endl;
  std::cout << "Testing VectorImageToImageAdaptor to extract a component from the vector image" << std::endl;


  const unsigned int componentToExtract = 2 * (Dimension -1);
  typedef itk::VectorImageToImageAdaptor< PixelType, Dimension > AdaptorType;
  AdaptorType::Pointer vectorImageToImageAdaptor = AdaptorType::New();
  vectorImageToImageAdaptor->SetExtractComponentIndex( componentToExtract );
  if( vectorImageToImageAdaptor->GetExtractComponentIndex() != componentToExtract )
    {
    std::cerr << "[FAILED]" << std::endl;
    }

  vectorImageToImageAdaptor->SetImage( vectorImage );
  vectorImageToImageAdaptor->Update();

  typedef itk::Image< PixelType , Dimension > AdaptedImageType;
  AdaptedImageType::IndexType index;
  index.Fill(10);

  if(   (vectorImageToImageAdaptor->GetPixel(index) !=  vectorImage->GetPixel( index )[componentToExtract])
     || (vectorImage->GetPixel( index )[componentToExtract] != componentToExtract ))
    {
    std::cerr << "[FAILED]" << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "[PASSED]" << std::endl;
    }
  }

  // Test with Region and Linear iterators...
  {
    // Create a  small image
    VectorImageType::Pointer vectorImage = VectorImageType::New();
    VectorImageType::IndexType start;
    itk::VariableLengthVector< PixelType > f( VectorLength );
    itk::VariableLengthVector< PixelType > ZeroPixel( VectorLength );
    ZeroPixel.Fill( itk::NumericTraits< PixelType >::Zero );
    for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i; }
    start.Fill(0);
    VectorImageType::SizeType  size;
    size.Fill(11);
    size[Dimension-1] = 5;
    unsigned long midCtr = 1;
    for (unsigned int i = 0; i < Dimension; i++) { midCtr *= size[i]; }
    VectorImageType::RegionType region( start, size );
    vectorImage->SetVectorLength( VectorLength );
    vectorImage->SetRegions( region );
    vectorImage->Allocate();
    vectorImage->FillBuffer( ZeroPixel );

    start.Fill(3);
    start[Dimension-1] = 2;
    size.Fill(4);
    size[Dimension-1] = 2;
    VectorImageType::RegionType subRegion( start, size );
    typedef itk::ImageRegionIterator< VectorImageType > ImageRegionIteratorType;
    ImageRegionIteratorType rit( vectorImage, subRegion );
    rit.GoToBegin();

    while( !rit.IsAtEnd() )
      {
      rit.Set( f );
      ++rit;
      }

    typedef itk::ImageRegionConstIterator< VectorImageType > ConstIteratorType;
    ConstIteratorType cit( vectorImage, vectorImage->GetBufferedRegion() );
    unsigned long ctr = 0;
    cit.Begin();
    midCtr /= 2;
    while( !cit.IsAtEnd() )
      {
      itk::VariableLengthVector< PixelType > value = cit.Get();
      ++cit;
      if( ctr == midCtr )
        {
        if( value != f )
          {
          std::cerr <<
            "ImageRegionConstIteratorTest on VectorImage [FAILED]" << std::endl;
          failed = true;
          }
        }
      ++ctr;
      }
    std::cout << "ImageRegionConstIteratorTest on VectorImage [PASSED]" << std::endl;


    {
    // Test itkImageLinearIteratorWithIndex
    typedef itk::ImageLinearConstIteratorWithIndex< VectorImageType > LinearConstIteratorType;
    typedef itk::ImageLinearIteratorWithIndex< VectorImageType > LinearIteratorType;

    LinearConstIteratorType lcit( vectorImage, vectorImage->GetBufferedRegion() );
    lcit.SetDirection( Dimension-1 );
    lcit.GoToBegin();
    itk::VariableLengthVector< PixelType > value;
    while( !lcit.IsAtEnd() )
      {
      while( !lcit.IsAtEndOfLine() )
        {
        value = lcit.Get();
        if( subRegion.IsInside( lcit.GetIndex() ) )
          {
          if( value!=f )
            {
            std::cerr <<
              "ImageLinearConstIteratorWithIndex on VectorImage [FAILED]" << std::endl;
            failed = true;
            }
          }
        else
          {
          if( value!=ZeroPixel )
            {
            std::cerr <<
              "ImageLinearConstIteratorWithIndex on VectorImage [FAILED]" << std::endl;
            failed = true;
            }
          }
        ++lcit;
        }
      lcit.NextLine();
      }

    VectorImageType::IndexType idx;
    idx.Fill(1);
    LinearIteratorType lit( vectorImage, vectorImage->GetBufferedRegion() );
    lit.SetIndex( idx );
    lit.Set( f );

    lcit.SetIndex( idx );
    value = lcit.Get();
    if( value != f )
      {
      std::cerr <<
        "ImageLinearConstIteratorWithIndex on VectorImage [FAILED]" << std::endl;
      failed = true;
      }

    std::cout << "ImageLinearConstIteratorWithIndex on VectorImage [PASSED]" << std::endl;
    std::cout << "ImageLinearIteratorWithIndex on VectorImage [PASSED]" << std::endl;
    }

  }

  }


  // Test IO support.
  {
  // Create an image using itk::Vector
  typedef itk::Vector< PixelType, VectorLength > VectorPixelType;
  typedef itk::Image< itk::Vector< PixelType, VectorLength >,
                                      Dimension > VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();
  VectorImageType::IndexType start;
  start.Fill(0);
  VectorImageType::SizeType size;
  size.Fill(5);
  VectorImageType::RegionType region( start, size );
  image->SetRegions( region );
  image->Allocate();

  typedef itk::ImageRegionIteratorWithIndex< VectorImageType > IteratorType;
  IteratorType it( image, region );
  it.GoToBegin();

  while( !it.IsAtEnd() )
    {
    VectorPixelType f;
    for (unsigned int i = 0; i < Dimension; i++)
      {
      f[i] = it.GetIndex()[i];
      f[Dimension+i] = it.GetIndex()[i];
      }
    it.Set( f );
    ++it;
    }

  typedef itk::ImageFileWriter< VectorImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( argv[2] );
  writer->Update();

  writer->SetFileName( argv[1] );
  writer->Update();
  }

  {
  // Now read it as a itk::VectorImage.
  typedef itk::VectorImage< PixelType, Dimension > VectorImageType;
  typedef itk::ImageFileReader< VectorImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  VectorImageType::Pointer vectorImage = reader->GetOutput();

  typedef itk::ImageRegionConstIteratorWithIndex< VectorImageType > IteratorType;
  IteratorType cit( vectorImage, vectorImage->GetBufferedRegion() );
  cit.GoToBegin();

  bool failed1 = false;
  while( !cit.IsAtEnd() )
    {
    for (unsigned int i = 0; i < Dimension; i++)
      {
      if (cit.Get()[i] != cit.GetIndex()[i] ||
          cit.Get()[i+Dimension] != cit.GetIndex()[i])
        {
        failed1 = true;
        }
      }
    ++cit;
    }

  if( failed1 )
    {
    std::cerr << "Read VectorImage [FAILED]" << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "Read VectorImage [PASSED]" << std::endl;
    }


  // Now write this out this VectorImage and read it again
  typedef itk::ImageFileWriter< VectorImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( vectorImage );
  writer->SetFileName(argv[1]);
  writer->Update();
  }


  {
  // Now read it as a itk::VectorImage.
  typedef itk::VectorImage< PixelType, Dimension > VectorImageType;
  typedef itk::ImageFileReader< VectorImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  VectorImageType::Pointer vectorImage = reader->GetOutput();

  typedef itk::ImageRegionConstIteratorWithIndex< VectorImageType > IteratorType;
  IteratorType cit( vectorImage, vectorImage->GetBufferedRegion() );
  cit.GoToBegin();

  bool failed1 = false;
  while( !cit.IsAtEnd() )
    {
    for (unsigned int i = 0; i < Dimension; i++)
      {
      if (cit.Get()[i] != cit.GetIndex()[i] ||
          cit.Get()[i+Dimension] != cit.GetIndex()[i])
        {
        failed1 = true;
        }
      }
    ++cit;
    }

  if( failed1 )
    {
    std::cerr << "Write VectorImage [FAILED]" << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "Write VectorImage [PASSED]" << std::endl;
    }


    {
    // Check support for Neighborhood Iterators
    //
    // 1. Test ConstNeighborhoodIterator
    //
    std::cout << "Testing ConstNeighborhoodIterator...." << std::endl;

    typedef itk::ConstNeighborhoodIterator< VectorImageType >
                                         ConstNeighborhoodIteratorType;
    ConstNeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);

    ConstNeighborhoodIteratorType::RegionType region
                            = vectorImage->GetBufferedRegion();
    ConstNeighborhoodIteratorType::SizeType size;
    size.Fill(4);
    ConstNeighborhoodIteratorType::IndexType index;
    index.Fill(1);
    region.SetIndex(index);
    region.SetSize( size );

    ConstNeighborhoodIteratorType cNit(radius, vectorImage, region);

    // Move Iterator to a point and see if it reads out the right value
    //
    unsigned int centerIndex = 1;
    for (unsigned int i = 0; i < Dimension; i++)
      { centerIndex *= (radius[i]*2+1); }
    centerIndex /= 2;

    ConstNeighborhoodIteratorType::IndexType location;
    for (unsigned int i = 0; i < Dimension; i++)
      {
      location[i] = i+1;
      }
    cNit.SetLocation( location );

    if( cNit.GetPixel(centerIndex) != vectorImage->GetPixel( location ) )
      {
      std::cerr << "  SetLocation [FAILED]" << std::endl;
      failed=true;
      }

    // Test GoToBegin()
    cNit.GoToBegin();
    if( cNit.GetPixel(centerIndex) != vectorImage->GetPixel( index ) )
      {
      std::cerr << "  GoToBegin [FAILED]" << std::endl;
      failed=true;
      }

    // Test GoToEnd()
    cNit.GoToEnd();
    --cNit;
    ConstNeighborhoodIteratorType::IndexType endIndex;
    endIndex.Fill(4);
    if( cNit.GetPixel(centerIndex) != vectorImage->GetPixel( endIndex ) )
      {
      std::cerr << "  GoToEnd [FAILED]" << std::endl;
      failed=true;
      }

    // Test IsAtEnd()
    if( !((++cNit).IsAtEnd()) )
      {
      std::cerr << "  IsAtEnd() [FAILED]" << std::endl;
      failed = true;
      }
    cNit.GoToBegin();
    unsigned int numPixelsTraversed = 1;
    for (unsigned int i = 0 ; i < Dimension; i++)
      { numPixelsTraversed *= size[i]; }
    while (! cNit.IsAtEnd())
      {
      ++cNit;
      --numPixelsTraversed;
      }
    if( numPixelsTraversed ) { std::cerr << "  IsAtEnd() [FAILED]" << std::endl; }

    // Test operator-
    --cNit;
    ConstNeighborhoodIteratorType::OffsetType offset;
    offset.Fill(1);
    cNit -= offset;
    itk::VariableLengthVector< PixelType > pixel = cNit.GetCenterPixel();
    itk::VariableLengthVector< PixelType > correctAnswer( VectorLength );
    correctAnswer.Fill( 3 );
    if( pixel != correctAnswer )
      {
      std::cerr << "  operator- [FAILED]" << std::endl;
      failed = true;
      }

    // Test GetNeighborhood()
    cNit.SetLocation( location );
    ConstNeighborhoodIteratorType::NeighborhoodType
                    neighborhood = cNit.GetNeighborhood();
    //const unsigned int neighborhoodSize = neighborhood.Size();
    //for( unsigned int i=0; i< neighborhoodSize; i++)
    //  { std::cout << neighborhood[i] << std::endl; }
    if( (neighborhood[0][0] != 0) || (neighborhood[0][2*Dimension-1] != (Dimension-1)))
      {
      std::cerr << "  GetNeighborhood() on ConstNeighborhoodIterator [FAILED]" << std::endl;
      failed = true;
      }



    //
    // 2. Test NeighborhoodIterator on VectorImage
    //

    std::cout << "Testing NeighborhoodIterator..." << std::endl;

    typedef itk::NeighborhoodIterator< VectorImageType > NeighborhoodIteratorType;
    NeighborhoodIteratorType nit(radius, vectorImage, region);
    nit.SetLocation( location );
    itk::VariableLengthVector< PixelType > p( VectorLength );
    p.Fill( 100.0 );
    nit.SetNext( 1, 1, p );

    // Test SetNext()
    NeighborhoodIteratorType::IndexType index1;
    index1 = location; index1[1] = location[1] + 1;
    nit.SetLocation( index1 );

    if( nit.GetCenterPixel() != p )
      {
      std::cerr << "  SetNext() [FAILED]" << std::endl;
      failed = true;
      }
    for (unsigned i = 0; i < Dimension; i++)
      {
      p[i] = p[Dimension + i] = (float)index1[i];
      }
    nit.SetCenterPixel( p );
    if( nit.GetCenterPixel() != p )
      {
      std::cerr << "  SetCenterPixel() [FAILED]" << std::endl;
      failed = true;
      }

    // Test SetNeighborhood() and GetPrevious()
    nit.SetLocation( index1 );
    nit.SetNeighborhood( neighborhood );
    for (unsigned i = 0; i < Dimension; i++)
      {
      p[i] = p[Dimension + i] = i + 1;
      }
    p[Dimension-1] = p[2*Dimension - 1] = Dimension - 1;
    if( nit.GetPrevious( Dimension-1, 1 ) != p )
      {
      std::cerr << "  SetNeighborhood() or GetPrevious() [FAILED]" << std::endl;
      failed = true;
      }

    if (Dimension == 3)
      {

      //
      // 3. Testing ConstShapedNeighborhoodIterator on VectorImage
      //

      // Go back to original image, where pixel values tell us the indices
      std::cout << "Testing ConstShapedNeighborhoodIterator on VectorImage..."
                                                                    << std::endl;
      reader->SetFileName( "dummy.nrrd");
      reader->SetFileName( argv[1] );
      reader->Update();
      vectorImage = reader->GetOutput();

      typedef itk::ConstShapedNeighborhoodIterator< VectorImageType >
                                     ConstShapedNeighborhoodIteratorType;
      ConstShapedNeighborhoodIteratorType cSnit( radius, vectorImage, region );
      cSnit.SetLocation( location );
      ConstShapedNeighborhoodIteratorType::OffsetType offset1;
      offset1[0] = 0; offset1[1] = 0; offset1[2] = 0;
      cSnit.ActivateOffset( offset1 ); //activate the center
      // activate the top plane
      offset1[0] = 0; offset1[1] = 0; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = 0; offset1[1] = 1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = 0; offset1[1] = -1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = 0; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = 1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = -1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = 0; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = 1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = -1; offset1[2] = -1;
      cSnit.ActivateOffset( offset1 );

      ConstShapedNeighborhoodIteratorType::IndexListType l
                                    = cSnit.GetActiveIndexList();
      ConstShapedNeighborhoodIteratorType::IndexListType::const_iterator
                                                          ali = l.begin();
      while (ali != l.end())
        {
        std::cout << *ali << " ";
        ++ali;
        }
      std::cout << std::endl;

      ConstShapedNeighborhoodIteratorType::ConstIterator ci = cSnit.Begin();
      while (! ci.IsAtEnd())
        {
        ConstShapedNeighborhoodIteratorType::OffsetType offset2 = ci.GetNeighborhoodOffset();
        if( (offset2[0] == -1) && (offset2[1]== -1) && (offset2[2]== -1) )
          {
          if( ci.GetNeighborhoodIndex() != 0 )
            {
            failed = true;
            std::cerr << "GetNeighborhoodOffset() on ConstShapedNeighborhoodIterato [FAILED]"
                                                                                << std::endl;
            }
          if( (ci.Get()[0]!=0) || (ci.Get()[1]!=1) || (ci.Get()[2]!=2) )
            {
            failed=true;
            std::cerr
              << "ConstShapedNeighborhoodIterator returned incorrect index [FAILED]"
                                                                          << std::endl;
            }
          }
        ci++;
        }

      //
      // 4. Test ShapedNeighborhoodIterator
      //
      typedef itk::ShapedNeighborhoodIterator< VectorImageType >
                                        ShapedNeighborhoodIteratorType;
      ShapedNeighborhoodIteratorType sNit( radius, vectorImage, region );

      offset1[0] = 0; offset1[1] = 0; offset1[2] = 0;
      sNit.ActivateOffset( offset1 ); //activate the center
      // activate the top plane
      offset1[0] = 0; offset1[1] = 0; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = 0; offset1[1] = 1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = 0; offset1[1] = -1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = 0; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = 1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = 1; offset1[1] = -1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = 0; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = 1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );
      offset1[0] = -1; offset1[1] = -1; offset1[2] = -1;
      sNit.ActivateOffset( offset1 );

      sNit.SetLocation( location );
      ShapedNeighborhoodIteratorType::Iterator shit = sNit.Begin();
      shit = sNit.Begin();
      p[0] = p[3] = 10; p[1] = p[4] = 20; p[2] = p[5] = 30;
      shit.Set( p );
      index[0]=location[0]-1; index[1]=location[1]-1; index[2]=location[2]-1;
      cNit.SetLocation( index );
      if( cNit.GetCenterPixel() != p )
        {
        std::cerr << "ShapedNeighborhoodIterator Set() [FAILED]" << std::endl;
        failed=true;
        }
      }


    }  // End Testing Neighborhood Iterators on VectorImage

  }

  if( failed )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
