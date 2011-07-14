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

#include "itkImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkQuasiNewtonObjectOptimizer.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethodImageSource.h"

/**
 *  This program tests the automatic estimation of parameter scales in optimizers
 *
 *
 */

int itkQuasiNewtonObjectOptimizerTest_Func( int argc,
                                           char* argv[],
                                           bool subtractMean )
{
  bool pass = true;

  const unsigned int dimension = 2;

  // Fixed Image Type
  typedef itk::Image<float,dimension>               FixedImageType;

  // Moving Image Type
  typedef itk::Image<float,dimension>               MovingImageType;

  // Size Type
  typedef MovingImageType::SizeType                 SizeType;


  // ImageSource
  typedef itk::testhelper::ImageRegistrationMethodImageSource<
                                  FixedImageType::PixelType,
                                  MovingImageType::PixelType,
                                  dimension >         ImageSourceType;

  // Transform Type
  typedef itk::AffineTransform< double, dimension > TransformType;
  typedef TransformType::ParametersType             ParametersType;

  // Optimizer Type
  typedef itk::QuasiNewtonObjectOptimizer                 OptimizerType;
  //typedef itk::AutomaticGradientDescentOptimizer    OptimizerType;

  // Metric Type
  typedef itk::NormalizedCorrelationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >    MetricType;

  // Interpolation technique
  typedef itk:: LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double >             InterpolatorType;

  // Registration Method
  typedef itk::ImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  typedef itk::CommandIterationUpdate<
                                    OptimizerType >      CommandIterationType;


  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  TransformType::Pointer      trasform      = TransformType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  unsigned long   numberOfIterations =   100;

  if( argc > 1 )
    {
    numberOfIterations = atol( argv[1] );
    std::cout << "numberOfIterations = " << numberOfIterations << std::endl;
    }

  FixedImageType::ConstPointer    fixedImage;
  MovingImageType::ConstPointer   movingImage;
  ImageSourceType::Pointer        imageSource;

  // Use image files as data for manual testing
  if( argc > 3 )
    {
    typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
    typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

    FixedImageReaderType::Pointer fixedImageReader   = FixedImageReaderType::New();
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

    fixedImageReader->SetFileName( argv[2] );
    movingImageReader->SetFileName( argv[3] );
    fixedImageReader->Update();
    movingImageReader->Update();

    fixedImage  = fixedImageReader->GetOutput();
    movingImage = movingImageReader->GetOutput();

    }
  else // Use imageSource to generate testing data
    {
    imageSource   = ImageSourceType::New();

    SizeType size;
    size[0] = 100;
    size[1] = 100;

    imageSource->GenerateImages( size );

    fixedImage    = imageSource->GetFixedImage();
    movingImage   = imageSource->GetMovingImage();
    }

  //
  // Connect all the components required for Registratio
  //
  registration->SetMetric(        metric        );
  //registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetFixedImage(    fixedImage    );
  registration->SetMovingImage(   movingImage   );
  registration->SetInterpolator(  interpolator  );

  // Select the Region of Interest over which the Metric will be computed
  // Registration time will be proportional to the number of pixels in this region.
  metric->SetFixedImageRegion( fixedImage->GetBufferedRegion() );

  // Turn on/off subtract mean flag.
  metric->SetSubtractMean( subtractMean );

  // Instantiate an Observer to report the progress of the Optimization
  CommandIterationType::Pointer iterationCommand = CommandIterationType::New();
  iterationCommand->SetOptimizer(  optimizer.GetPointer() );

  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetMaximize(false);

  // Start from an Identity transform (in a normal case, the user
  // can probably provide a better guess than the identity...
  transform->SetIdentity();
  registration->SetInitialTransformParameters( transform->GetParameters() );

  /*// Testing optimizer parameter estimator
  typedef itk::OptimizerParameterEstimator< FixedImageType,
                                        MovingImageType,
                                        TransformType > OptimizerParameterEstimatorType;
  OptimizerParameterEstimatorType::Pointer parameterEstimator = OptimizerParameterEstimatorType::New();

  parameterEstimator->SetFixedImage(fixedImage);
  parameterEstimator->SetMovingImage(movingImage);
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromShift);
  //parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromJacobian);

  //parameterEstimator->DebugOn();
  //optimizer->DebugOn();

  optimizer->SetOptimizerHelper( parameterEstimator );
  // Setting optimizer parameter done
  */
  // Initialize the internal connections of the registration method.
  // This can potentially throw an exception
  try
    {
    registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    pass = false;
    }

  if( argc > 3 ) // Use image files as testing data
    {
    ParametersType finalParameters  = registration->GetLastTransformParameters();
    std::cout << "finalParameters = " << finalParameters << std::endl;
    std::cout << "Test PASSED." << std::endl;
    return EXIT_SUCCESS;
    }

  ParametersType actualParameters = imageSource->GetActualParameters();
  ParametersType finalParameters  = registration->GetLastTransformParameters();

  const unsigned int numbeOfParameters = actualParameters.Size();

  // We know that for the Affine transform the Translation parameters are at
  // the end of the list of parameters.
  const unsigned int offsetOrder = finalParameters.Size()-actualParameters.Size();

  const double tolerance = 1.0;  // equivalent to 1 pixel.
  std::cout << "Estimated scales = " << optimizer->GetScales() << std::endl;
  std::cout << "finalParameters = " << finalParameters << std::endl;
  std::cout << "actualParameters = " << actualParameters << std::endl;

  for(unsigned int i=0; i<numbeOfParameters; i++)
    {
    // the parameters are negated in order to get the inverse transformation.
    // this only works for comparing translation parameters....
    std::cout << finalParameters[i+offsetOrder] << " == " << -actualParameters[i] << std::endl;
    if( vnl_math_abs ( finalParameters[i+offsetOrder] - (-actualParameters[i]) ) > tolerance )
      {
      std::cout << "Tolerance exceeded at component " << i << std::endl;
      pass = false;
      }
    }

  if( !pass )
    {
    std::cout << "Test FAILED." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test PASSED." << std::endl;
  return EXIT_SUCCESS;


}

int itkQuasiNewtonObjectOptimizerTest( int argc, char* argv[] )
{
  // test metric without factoring out the mean.
  int fail1 = itkQuasiNewtonObjectOptimizerTest_Func( argc, argv, false );

  // test metric with factoring out the mean.
  int fail2 = itkQuasiNewtonObjectOptimizerTest_Func( argc, argv, true );

  if( fail1 || fail2 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}
