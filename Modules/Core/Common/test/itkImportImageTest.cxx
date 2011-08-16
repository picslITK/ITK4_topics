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
#include "itkImageRegionIterator.h"
#include "itkShrinkImageFilter.h"
#include "itkImportImageFilter.h"

int itkImportImageTest(int, char* [] )
{
  // Create a C-array to hold an image
  short *rawImage = new short[8*12];
  for (int i=0; i < 8*12; i++)
    {
    rawImage[i] = i;
    }

  // typdefs to simplify the syntax
  typedef itk::ImportImageFilter<short, 2>          ImportImageFilter;
  typedef itk::Image<short, 2>   ShortImage;

  // Create an ImportImageFilter filter
  ImportImageFilter::Pointer import;
  import = ImportImageFilter::New();

  itk::ImageRegion<2>         region;
  itk::ImageRegion<2>::IndexType  index = {{0, 0}};
  itk::ImageRegion<2>::SizeType   size = {{8, 12}};

  region.SetSize( size );
  region.SetIndex( index );

  import->SetRegion( region );
  import->SetImportPointer( rawImage, 8*12, true);

  // Create another filter
  itk::ShrinkImageFilter<ImportImageFilter::OutputImageType, ShortImage >::Pointer shrink;
  shrink = itk::ShrinkImageFilter<ImportImageFilter::OutputImageType, ShortImage>::New();
  shrink->SetInput( import->GetOutput() );
  shrink->SetShrinkFactors(2); //Also tested with factors 3 and 4, with 12x12 image
  try
    {
    shrink->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception detected: "  << e.GetDescription();
    return EXIT_FAILURE;
    }



  // Test the SetVectorMacros and GetVectorMacros
  const double data[2] = {1.0,1.0};
  import->SetSpacing(data);

  const float data2[2] = {1.0,1.0};
  import->SetSpacing(data2);

  const double * spacingValue = import->GetSpacing().GetDataPointer();
  std::cout << "import->GetSpacing(): " << spacingValue << std::endl;

  const double data3[2] = {1.0,1.0};
  import->SetOrigin(data3);

  const float data4[2] = {1.0,1.0};
  import->SetOrigin(data4);

  const double * originValue = import->GetOrigin().GetDataPointer();
  std::cout << "import->GetOrigin(): " << originValue << std::endl;

  //
  // The rest of this code determines whether the shrink code produced
  // the image we expected.
  //
  ShortImage::RegionType requestedRegion;
  requestedRegion = shrink->GetOutput()->GetRequestedRegion();

  itk::ImageRegionIterator<ShortImage>
    iterator2(shrink->GetOutput(), requestedRegion);

  bool passed = true;
  for (; !iterator2.IsAtEnd(); ++iterator2)
    {
    std::cout << "Pixel " << iterator2.GetIndex() << " = " << iterator2.Get() << std::endl;
    if ( iterator2.Get() !=
         itk::Math::RoundHalfIntegerUp<short>(
            static_cast<float>( (shrink->GetShrinkFactors()[0] * iterator2.GetIndex()[0] +
                                 shrink->GetShrinkFactors()[0]/2) +
            (region.GetSize()[0] * ((shrink->GetShrinkFactors()[1]/2) +
                                    (shrink->GetShrinkFactors()[0] * iterator2.GetIndex()[1]))))))
      {
         std::cout << " iterator2.GetIndex() Get() " << iterator2.GetIndex() <<  " " << iterator2.Get()
                   << " compare value " << itk::Math::RoundHalfIntegerUp<short>(
            static_cast<float>( (shrink->GetShrinkFactors()[0] * iterator2.GetIndex()[0] +
                                 shrink->GetShrinkFactors()[0]/2) +
            (region.GetSize()[0] * ((shrink->GetShrinkFactors()[1]/2) +
                                    (shrink->GetShrinkFactors()[0] * iterator2.GetIndex()[1]))))) << "\n";
      passed = false;
      }
    }

  if (passed)
    {
    std::cout << "ImportImageFilter test passed." << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "ImportImageFilter test failed." << std::endl;
    return EXIT_FAILURE;
    }


}
