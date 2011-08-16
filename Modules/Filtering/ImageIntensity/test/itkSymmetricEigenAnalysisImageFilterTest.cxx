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
#include "itkSymmetricSecondRankTensor.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"


int itkSymmetricEigenAnalysisImageFilterTest(int, char* [] )
{

  // Define the dimension of the images
  const unsigned int myDimension = 3;

  // Define the symmetric tensor pixel type
  typedef itk::SymmetricSecondRankTensor< float, myDimension > myTensorType;

  // Declare the types of the images
  typedef itk::Image< myTensorType, myDimension >  myInputImageType;

  // Define the type for storing the eigen-value
  typedef itk::FixedArray< float, myDimension >  myValueArray;

  // Declare the types of the output images
  typedef itk::Image< myValueArray, myDimension >  myOutputImageType;

  // Declare the type of the index to access images
  typedef itk::Index<myDimension>             myIndexType;

  // Declare the type of the size
  typedef itk::Size<myDimension>              mySizeType;

  // Declare the type of the Region
  typedef itk::ImageRegion<myDimension>        myRegionType;

  // Create the image
  myInputImageType::Pointer inputImage  = myInputImageType::New();


  // Define their size, and start index
  mySizeType size;
  size[0] = 8;
  size[1] = 8;
  size[2] = 8;

  myIndexType start;
  start.Fill(0);

  myRegionType region;
  region.SetIndex( start );
  region.SetSize( size );

  // Initialize Image A
  inputImage->SetLargestPossibleRegion( region );
  inputImage->SetBufferedRegion( region );
  inputImage->SetRequestedRegion( region );
  inputImage->Allocate();

  // Declare Iterator type for the input image
  typedef itk::ImageRegionIteratorWithIndex<
                     myInputImageType>  myIteratorType;

  // Create one iterator for the Input Image A (this is a light object)
  myIteratorType it( inputImage, inputImage->GetRequestedRegion() );

  myTensorType tensorValue;

  tensorValue(0,0) = 19.0;
  tensorValue(0,1) = 23.0;
  tensorValue(0,2) = 29.0;
  tensorValue(1,1) = 31.0;
  tensorValue(1,2) = 37.0;
  tensorValue(2,2) = 39.0;

  it.GoToBegin();

  // Initialize the content of Image A
  while( !it.IsAtEnd() )
    {
    it.Set( tensorValue );
    ++it;
    }


  // Declare the type for the filter
  typedef itk::SymmetricEigenAnalysisImageFilter<
                                     myInputImageType,
                                     myOutputImageType
                                               >  myFilterType;



  // Create a  Filter
  myFilterType::Pointer filter = myFilterType::New();
  filter->SetDimension( myTensorType::Dimension );

  // Connect the input images
  filter->SetInput( inputImage );


  // Execute the filter
  filter->Update();
  filter->SetFunctor(filter->GetFunctor());

  // Get the Smart Pointer to the Filter Output
  // It is important to do it AFTER the filter is Updated
  // Because the object connected to the output may be changed
  // by another during GenerateData() call
  myOutputImageType::Pointer outputImage = filter->GetOutput();

  // Declare Iterator type for the output image
  typedef itk::ImageRegionIteratorWithIndex<
                                 myOutputImageType>  myOutputIteratorType;

  // Create an iterator for going through the output image
  myOutputIteratorType itg( outputImage,
                            outputImage->GetRequestedRegion() );

  //  Print the content of the result image
  std::cout << " Result " << std::endl;
  itg.GoToBegin();
  while( !itg.IsAtEnd() )
    {
    std::cout << itg.Get();
    ++itg;
    }


  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;

}




