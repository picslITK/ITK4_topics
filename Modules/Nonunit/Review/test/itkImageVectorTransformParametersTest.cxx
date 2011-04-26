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

#include "itkImageVectorTransformParameters.h"

using namespace itk;

namespace{

typedef double                                        ValueType;
const   unsigned int                                  ImageDimension = 2;
const   unsigned int                                  VectorDimension = 4;
typedef itk::Vector< ValueType, VectorDimension >     VectorPixelType;
typedef itk::Image< VectorPixelType, ImageDimension > ImageVectorType;
typedef ImageVectorType::Pointer                      ImageVectorPointer;
typedef ImageVectorType::RegionType                   RegionType;
typedef RegionType::SizeType                          SizeType;
typedef ImageVectorType::IndexType                    IndexType;
typedef ImageVectorType::PixelContainer               VectorPixelContainer;
typedef ImageVectorTransformParameters< ValueType,
                                        VectorDimension,
                                        ImageDimension >
                                          ImageVectorTransformParametersType;
}

int testMemoryAccess( ImageVectorTransformParametersType& imageVectorParams,
                      ImageVectorPointer imageVector,
                      int dimLength )
{
  int result = EXIT_SUCCESS;

  for (int y = 0; y < dimLength; y++)
    {
    for (int x = 0; x < dimLength; x++)
      {
      IndexType index;
      index[0] = x;
      index[1] = y;

      // The image index returns a N-dim vector, so have to check each
      // element against the values returned by parameter object.
      unsigned long offset = (x + y * dimLength) * VectorDimension;
      VectorPixelType vectorpixel = imageVector->GetPixel( index );
      for(unsigned int ind=0; ind < VectorDimension; ind++)
        {
        ValueType paramsValue = imageVectorParams[offset+ind];
        if( vectorpixel[ind] != paramsValue )
          {
          std::cout << "VectorImage pixel value does not match params value."
                    << "vectorpixel[" << ind << "]: " << vectorpixel[ind]
                    << std::endl
                    << "imageVectorParams[" << offset+ind << "]: "
                    << paramsValue << std::endl;
          result = EXIT_FAILURE;
          }
        }
      }
    }
  return result;
}

/******************************************************/

int itkImageVectorTransformParametersTest(int, char *[])
{
  int result = EXIT_SUCCESS;

  ImageVectorPointer imageVector = ImageVectorType::New();

  IndexType start;
  start.Fill( 0 );

  SizeType size;
  const int dimLength = 3;
  size.Fill( dimLength );

  RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  imageVector->SetRegions( region );
  imageVector->Allocate();

  ImageVectorType::PointType     origin;
  ImageVectorType::SpacingType   spacing;

  origin.Fill( 0.0 );
  spacing.Fill( 1.0 );

  imageVector->SetOrigin( origin );
  imageVector->SetSpacing( spacing );

  ValueType vectorinitvalues[VectorDimension] = {0.0, 0.1, 0.2, 0.3};
  VectorPixelType vectorvalues(vectorinitvalues);

  //
  // Fill up the image values with the function
  //
  //   Intensity = f(x,y) = x + 3 * y
  //
  //
  for (int y = 0; y < dimLength; y++)
    {
    for (int x = 0; x < dimLength; x++)
      {
      IndexType index;
      index[0] = x;
      index[1] = y;

      const ValueType value = x + y * dimLength;

      VectorPixelType & vectorpixel = imageVector->GetPixel( index );
      vectorpixel.Fill( value );
      vectorpixel += vectorvalues;

      std::cout << value << " ";
      }
    std::cout << std::endl;
    }

  ImageVectorTransformParametersType imageVectorParams;
  imageVectorParams.SetParameterImage( imageVector );

  result = testMemoryAccess( imageVectorParams, imageVector, dimLength );

  return result;
}
