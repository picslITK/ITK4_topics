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
#include "itkCropImageFilter.h"
#include "itkFileOutputWindow.h"
#include "itkStreamingImageFilter.h"
#include "itkFilterWatcher.h"

int itkCropImageFilterTest(int, char* [] )
{
  itk::FileOutputWindow::Pointer fow = itk::FileOutputWindow::New();
  fow->SetInstance(fow);


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

  std::cout << "Input Image: " << if2 << std::endl;

  // Create a filter
  itk::CropImageFilter< ShortImage, ShortImage >::Pointer extract;
  extract = itk::CropImageFilter< ShortImage, ShortImage >::New();
  FilterWatcher watcher(extract);

  extract->SetInput( if2 );

  ShortImage::RegionType requestedRegion;

  ShortImage::SizeType   extractSize = {{8, 12}};
  extractSize[0] = 1; extractSize[1] = 1;
  extract->SetBoundaryCropSize(extractSize);
  extract->SetUpperBoundaryCropSize(extractSize);
  extract->SetLowerBoundaryCropSize(extractSize);
  extract->UpdateLargestPossibleRegion();
  requestedRegion = extract->GetOutput()->GetRequestedRegion();

  if (extract->GetOutput()->GetLargestPossibleRegion().GetSize()[0] != 6
      || extract->GetOutput()->GetLargestPossibleRegion().GetSize()[1] != 10)
    {
      return EXIT_FAILURE;
    }

  if (extract->GetOutput()->GetLargestPossibleRegion().GetIndex()[0] != 1
      || extract->GetOutput()->GetLargestPossibleRegion().GetIndex()[1] != 1)
    {
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
