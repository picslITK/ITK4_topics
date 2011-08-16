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
#include "itkImage.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkConstantPadImageFilter.h"
#include "itkFileOutputWindow.h"
#include "itkStreamingImageFilter.h"
#include "itkFilterWatcher.h"

int itkConstantPadImageTest(int, char* [] )
{
  itk::FileOutputWindow::Pointer fow = itk::FileOutputWindow::New();
  fow->SetInstance(fow);

  int nextVal;

  // typedefs to simplify the syntax
  typedef itk::Image<short, 2>   SimpleImage;
  SimpleImage::Pointer simpleImage = SimpleImage::New();
  std::cout << "Simple image spacing: " << simpleImage->GetSpacing()[0] << ", "
            << simpleImage->GetSpacing()[1] << std::endl;

  // typedefs to simplify the syntax
  typedef itk::Image<short, 2>   ShortImage;

  // Test the creation of an image with native type
  ShortImage::Pointer if2 = ShortImage::New();

  // fill in an image
  ShortImage::IndexType  index = {{0, 0}};
  ShortImage::SizeType   size = {{8, 12}};
  ShortImage::RegionType region;
  int row, column;
  region.SetSize( size );
  region.SetIndex( index );
  if2->SetLargestPossibleRegion( region );
  if2->SetBufferedRegion( region );
  if2->Allocate();

  itk::ImageRegionIterator<ShortImage> iterator(if2, region);

  short i=0;
  for (; !iterator.IsAtEnd(); ++iterator, ++i)
    {
    iterator.Set( i );
    }

  // Create a filter
  itk::ConstantPadImageFilter< ShortImage, ShortImage >::Pointer constantPad;
  constantPad = itk::ConstantPadImageFilter< ShortImage, ShortImage >::New();
  FilterWatcher watch(constantPad);
  constantPad->SetInput( if2 );

  typedef ShortImage::SizeValueType   SizeValueType;
  typedef ShortImage::IndexValueType  IndexValueType;

  SizeValueType upperfactors[2] = { 0, 0};
  SizeValueType lowerfactors[2] = { 0, 0};

  constantPad->SetConstant(13);
  // check the method using the SizeType rather than the simple table type.
  ShortImage::SizeType stfactors;
  stfactors.Fill( 0 );
  constantPad->SetPadLowerBound(stfactors);
  constantPad->SetPadUpperBound(stfactors);
  constantPad->SetPadBound(stfactors);
  constantPad->UpdateLargestPossibleRegion();

  std::cout << constantPad << std::endl;
  std::cout << "Input spacing: " << if2->GetSpacing()[0] << ", "
            << if2->GetSpacing()[1] << std::endl;
  std::cout << "Output spacing: " << constantPad->GetOutput()->GetSpacing()[0]
            << ", "
            << constantPad->GetOutput()->GetSpacing()[1] << std::endl;


  ShortImage::RegionType requestedRegion;
  bool passed;

  // CASE 1
  lowerfactors[0] = 1; lowerfactors[1] = 2;
  upperfactors[0] = 3; upperfactors[1] = 4;
  constantPad->SetPadLowerBound(lowerfactors);
  constantPad->SetPadUpperBound(upperfactors);
  constantPad->UpdateLargestPossibleRegion();
  requestedRegion = constantPad->GetOutput()->GetRequestedRegion();

  itk::ImageRegionIterator<ShortImage>
    iteratorIn1(constantPad->GetOutput(), requestedRegion);

  passed = true;
  size = requestedRegion.GetSize();
  index = requestedRegion.GetIndex();
  if ((index[0] != (0 - (IndexValueType) lowerfactors[0]))
      || (index[1] != (0 - (IndexValueType) lowerfactors[1]))
      || (size[0] != (8 + lowerfactors[0] + upperfactors[0]))
      || (size[1] != (12 + lowerfactors[1] + upperfactors[1])))
    {
    passed = false;
    }
  else
    {
    for (; !iteratorIn1.IsAtEnd(); ++iteratorIn1)
      {
      row = iteratorIn1.GetIndex()[0];
      column = iteratorIn1.GetIndex()[1];
      if ((row < 0) || (row>7) || (column < 0) || (column > 11))
        {
        if ( iteratorIn1.Get() != 13 )
          {
          passed = false;
          }
        }
      else
        {
        nextVal = 8*column+row;
        if (iteratorIn1.Get() != nextVal)
          {
          std::cout << "Error: (" << row << ", " << column
                    << "), expected " << nextVal << " got "
                    << iteratorIn1.Get() << std::endl;
          passed = false;
          }
        }
      }
    }

  if (passed)
    {
    std::cout << "constantPadImageFilter case 1 passed." << std::endl;
    }
  else
    {
    std::cout << "constantPadImageFilter case 1 failed." << std::endl;
    return EXIT_FAILURE;
    }


  // CASE 2
  lowerfactors[0] = 10;
  upperfactors[1] = 15;
  constantPad->SetPadLowerBound(lowerfactors);
  constantPad->SetPadUpperBound(upperfactors);

  // Create a stream
  itk::StreamingImageFilter< ShortImage, ShortImage >::Pointer stream;
  stream = itk::StreamingImageFilter< ShortImage, ShortImage >::New();
  stream->SetInput( constantPad->GetOutput() );
  stream->SetNumberOfStreamDivisions(2);


  if ((constantPad->GetPadUpperBound()[0] != upperfactors[0])
      || (constantPad->GetPadUpperBound()[1] != upperfactors[1])
      || (constantPad->GetPadLowerBound()[0] != lowerfactors[0])
      || (constantPad->GetPadLowerBound()[1] != lowerfactors[1]))
    {
    passed = false;
    }
  else
    {
    stream->UpdateLargestPossibleRegion();
    requestedRegion = stream->GetOutput()->GetRequestedRegion();

    itk::ImageRegionIterator<ShortImage>
      iteratorIn2(stream->GetOutput(), requestedRegion);

    passed = true;
    size = requestedRegion.GetSize();
    index = requestedRegion.GetIndex();
    if ((index[0] != (0 - (IndexValueType) lowerfactors[0]))
        || (index[1] != (0 - (IndexValueType) lowerfactors[1]))
        || (size[0] != (8 + lowerfactors[0] + upperfactors[0]))
        || (size[1] != (12 + lowerfactors[1] + upperfactors[1])))
      {
      passed = false;
      }
    else
      {
      for (; !iteratorIn2.IsAtEnd(); ++iteratorIn2)
        {
        row = iteratorIn2.GetIndex()[0];
        column = iteratorIn2.GetIndex()[1];
        if ((row < 0) || (row>7) || (column < 0) || (column > 11))
          {
          if ( iteratorIn2.Get() != 13 )
            {
            passed = false;
            }
          }
        else
          {
          nextVal = 8*column+row;
          if (iteratorIn2.Get() != nextVal)
            {
            std::cout << "Error: (" << row << ", " << column
                      << "), expected " << nextVal << " got "
                      << iteratorIn2.Get() << std::endl;
            passed = false;
            }
          }
        }
      }
    }

  if (passed)
    {
    std::cout << "constantPadImageFilter case 2 passed." << std::endl;
    }
  else
    {
    std::cout << "constantPadImageFilter case 2 failed." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
