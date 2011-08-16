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
#include "itkNotImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"


int itkNotImageFilterTest(int, char* [] )
{

  // Define the dimension of the images
  const unsigned int myDimension = 3;

  // Declare the types of the images
  typedef bool myPixelType;
  typedef itk::Image<myPixelType, myDimension>  myImageType1;
  typedef itk::Image<myPixelType, myDimension>  myImageType3;

  // Declare the type of the index to access images
  typedef itk::Index<myDimension>         myIndexType;

  // Declare the type of the size
  typedef itk::Size<myDimension>          mySizeType;

  // Declare the type of the Region
  typedef itk::ImageRegion<myDimension>        myRegionType;

  // Declare the type for the ADD filter
  typedef itk::NotImageFilter<
    myImageType1,
    myImageType3  >       myFilterType;

  // Declare the pointers to images
  typedef myImageType1::Pointer   myImageType1Pointer;
  typedef myImageType3::Pointer   myImageType3Pointer;
  typedef myFilterType::Pointer   myFilterTypePointer;

  // Create an image
  myImageType1Pointer inputImageA  = myImageType1::New();

  // Define their size, and start index
  mySizeType size;
  size[0] = 2;
  size[1] = 2;
  size[2] = 2;

  myIndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  myRegionType region;
  region.SetIndex( start );
  region.SetSize( size );

  // Initialize Image A
  inputImageA->SetLargestPossibleRegion( region );
  inputImageA->SetBufferedRegion( region );
  inputImageA->SetRequestedRegion( region );
  inputImageA->Allocate();

  // Declare Iterator types apropriated for each image
  typedef itk::ImageRegionIteratorWithIndex<myImageType1>  myIteratorType1;
  typedef itk::ImageRegionIteratorWithIndex<myImageType3>  myIteratorType3;

  // Create one iterator for Image A (this is a light object)
  myIteratorType1 it1( inputImageA, inputImageA->GetBufferedRegion() );
  it1.GoToBegin();

  // Initialize the content of Image A
  while( !it1.IsAtEnd() )
  {
    it1.Set( true );
    std::cout << static_cast<itk::NumericTraits<myPixelType>::PrintType>(it1.Get()) << std::endl;
    ++it1;
  }

  // Create a NOT Filter
  myFilterTypePointer filter = myFilterType::New();


  // Connect the input images
  filter->SetInput( inputImageA );

  // Get the Smart Pointer to the Filter Output
  myImageType3Pointer outputImage = filter->GetOutput();


  // Execute the filter
  filter->Update();
  filter->SetFunctor(filter->GetFunctor());

  // Create an iterator for going through the image output
  myIteratorType3 it3(outputImage, outputImage->GetBufferedRegion());
  it3.GoToBegin();

  //  Print the content of the result image
  std::cout << " Result " << std::endl;
  while( !it3.IsAtEnd() )
  {
    std::cout << static_cast<itk::NumericTraits<myPixelType>::PrintType>(it3.Get()) << std::endl;
    ++it3;
  }


  // All objects should be automatically destroyed at this point
  std::cout << "Test PASSED !" << std::endl;

  return EXIT_SUCCESS;

}




