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

#include "itkGaussianDerivativeImageFunction.h"
#include "itkImage.h"

int itkGaussianDerivativeImageFunctionTest(int, char* [] )
{
  const unsigned int Dimension = 2;
  typedef float  PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::GaussianDerivativeImageFunction< ImageType > DoGFunctionType;

  // Create and allocate the image
  ImageType::Pointer      image = ImageType::New();
  ImageType::SizeType     size;
  ImageType::IndexType    start;
  ImageType::RegionType   region;

  size[0] = 50;
  size[1] = 50;

  start.Fill( 0 );
  region.SetIndex( start );
  region.SetSize( size );

  image->SetRegions( region );
  image->Allocate();

  ImageType::PixelType initialValue = 0;
  image->FillBuffer( initialValue );

  // Fill the image with a straight line
  for(unsigned int i=0;i<50;i++)
  {
    ImageType::IndexType ind;
    ind[0]=i;
    ind[1]=25;
    image->SetPixel(ind,1);
    ind[1]=26;
    image->SetPixel(ind,1);
  }

  // Test the derivative of Gaussian image function
  DoGFunctionType::Pointer DoG = DoGFunctionType::New();
  DoG->SetInputImage( image );

  std::cout << "Testing Set/GetSigma(): ";

  DoG->SetSigma(2.0);
  const double* sigma = DoG->GetSigma();
  for(unsigned int i=0;i<Dimension;i++)
  {
    if( sigma[i] !=  2.0)
    {
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }
  }
  std::cout << "[PASSED] " << std::endl;


  std::cout << "Testing Set/GetExtent(): ";

  DoG->SetExtent(4.0);
  const double* ext = DoG->GetExtent();
  for(unsigned int i=0;i<Dimension;i++)
  {
    if( ext[i] !=  4.0)
    {
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }
  }
  std::cout << "[PASSED] " << std::endl;

  std::cout << "Testing consistency within Index/Point/ContinuousIndex: ";
  itk::Index<2>   index;
  index.Fill(25);
  DoGFunctionType::OutputType  gradient_index;
  gradient_index = DoG->EvaluateAtIndex( index );

  DoGFunctionType::PointType pt;
  pt[0]=25.0;
  pt[1]=25.0;
  DoGFunctionType::OutputType  gradient_point;
  gradient_point = DoG->Evaluate( pt );


  DoGFunctionType::ContinuousIndexType continuousIndex;
  continuousIndex.Fill(25);
  DoGFunctionType::OutputType  gradient_continuousIndex;
  gradient_continuousIndex = DoG->EvaluateAtContinuousIndex( continuousIndex );

  if( gradient_index !=  gradient_point
     || gradient_index != gradient_continuousIndex)
    {
    std::cerr << "[FAILED] : " << gradient_index << " : "
              << gradient_point << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "[PASSED] " << std::endl;
  gradient_point.Normalize(); // normalize the vector;

  std::cout << "Testing Evaluate() : ";

  if( (gradient_point[0] > 0.1)  ||
      (vcl_fabs(gradient_point[1]+1.0)> 10e-4)
    )
    {
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "[PASSED] " << std::endl;

  pt[0]=25.0;
  pt[1]=26.0;
  gradient_point = DoG->Evaluate( pt );

  gradient_point.Normalize(); // normalize the vector;

  std::cout << "Testing Evaluate() : ";

  if( (gradient_point[0] > 0.1)  ||
      (vcl_fabs(gradient_point[1]-1.0)> 10e-4)
    )
    {
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "[PASSED] " << std::endl;

  std::cout << "Testing Gaussian Derivative Spatial Function:";

  typedef itk::GaussianDerivativeSpatialFunction<double,1>  GaussianDerivativeFunctionType;
  GaussianDerivativeFunctionType::Pointer f = GaussianDerivativeFunctionType::New();

  f->SetScale(1.0);
  if(f->GetScale() != 1.0)
    {
    std::cerr << "Get Scale : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  f->SetNormalized(true);
  if(!f->GetNormalized())
    {
    std::cerr << "GetNormalized : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  GaussianDerivativeFunctionType::ArrayType s;
  s[0] = 1.0;
  f->SetSigma(s);
  if(f->GetSigma()[0] != 1.0)
    {
    std::cerr << "GetSigma : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  GaussianDerivativeFunctionType::ArrayType m;
  m[0] = 0.0;
  f->SetMean(m);
  if(f->GetMean()[0] != 0.0)
    {
    std::cerr << "GetMean : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  f->SetDirection(0);
  if(f->GetDirection() != 0)
    {
    std::cerr << "GetDirection : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  GaussianDerivativeFunctionType::InputType point;
  point[0] = 0.0;

  if(f->Evaluate(point) != 0.0)
    {
    std::cerr << "Evaluate: [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << f << std::endl;
  std::cout << "[PASSED] " << std::endl;

  std::cout << "GaussianDerivativeImageFunctionTest: [DONE] " << std::endl;
  return EXIT_SUCCESS;
}

