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

#include "itkAffineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

#include "itkAutomaticGradientDescentOptimizer.h"
#include "itkQuasiNewtonOptimizer.h"
#include "itkOptimizerParameterEstimator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"

#include "itkTextOutput.h"
#include "itkImageRegionIterator.h"
#include "itkCommandIterationUpdate.h"
#include "itkSimpleMultiResolutionImageRegistrationUI.h"

/**
 *  This program test the
 *  itk::QuasiNewtonOptimizer and
 *  itk::MultiResolutionImageRegistrationMethod class
 *
 *  This file tests the combination of:
 *   - NormalizedCorrelationImageToImageMetric
 *   - AffineTransform
 *   - QuasiNewtonOptimizer
 *   - OptimizerParameterEstimator
 *   - RecursiveMultiResolutionPyramidImageFilter
 *
 */
template< unsigned int ImageDimension >
int itkOptimizerParameterEstimatorTest_1_Dimension
    (int argc, char* argv[] );

int itkOptimizerParameterEstimatorTest_1(int argc, char* argv[] )
{
  int ret = EXIT_SUCCESS;
  if (argc < 5)
    {
    std::cout << "Usage: " << argv[0]
              << " ImageDimension FixedImage MovingImage TransformFile" << std::endl;
    return ret;
    }

  int dimension = atoi(argv[1]);

  if (dimension == 2)
    {
    ret = itkOptimizerParameterEstimatorTest_1_Dimension<2>(argc, argv);
    }
  else if (dimension == 3)
    {
    ret = itkOptimizerParameterEstimatorTest_1_Dimension<3>(argc, argv);
    }
  else
    {
    std::cout << "Only 2 or 3 dimensional images are supported." << std::endl;
    ret = EXIT_FAILURE;
    }
  return ret;

}
template< unsigned int ImageDimension >
int itkOptimizerParameterEstimatorTest_1_Dimension
    (int argc, char* argv[] )
{
  unsigned long   numberOfIterations =   200;

  itk::OutputWindow::SetInstance(itk::TextOutput::New().GetPointer());

  const unsigned int dimension = ImageDimension;

  typedef float  PixelType;

  // Fixed Image Type
  typedef itk::Image<PixelType,dimension>               FixedImageType;

  // Moving Image Type
  typedef itk::Image<PixelType,dimension>               MovingImageType;

  // Transform Type
  typedef itk::AffineTransform< double,dimension >  TransformType;

  // Optimizer Type
  typedef itk::QuasiNewtonOptimizer                     OptimizerType;

  // Metric Type
  typedef itk::NormalizedCorrelationImageToImageMetric<
                                    FixedImageType,
                                    MovingImageType >    MetricType;

  // Interpolation technique
  typedef itk:: LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double          >    InterpolatorType;

  // Fixed Image Pyramid Type
  typedef itk::RecursiveMultiResolutionPyramidImageFilter<
                                    FixedImageType,
                                    FixedImageType  >    FixedImagePyramidType;

  // Moving Image Pyramid Type
  typedef itk::RecursiveMultiResolutionPyramidImageFilter<
                                    MovingImageType,
                                    MovingImageType  >   MovingImagePyramidType;


  // Registration Method
  typedef itk::MultiResolutionImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  typedef itk::CommandIterationUpdate<
                                    OptimizerType >      CommandIterationType;

  /*********************************************************
   * Set up the two input images.
   * One image scaled and shifted with respect to the other.
   **********************************************************/
  typename FixedImageType::Pointer     fixedImage;
  typename MovingImageType::Pointer    movingImage;

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  typename FixedImageReaderType::Pointer fixedImageReader   = FixedImageReaderType::New();
  typename MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName( argv[2] );
  movingImageReader->SetFileName( argv[3] );
  fixedImageReader->Update();
  movingImageReader->Update();

  fixedImage  = fixedImageReader->GetOutput();
  movingImage = movingImageReader->GetOutput();

  typename RegistrationType::ScheduleType  fixedImageSchedule;
  typename RegistrationType::ScheduleType  movingImageSchedule;

  /* The registration run invokes SetNumberOfLevels to specify
   * the number of computation levels */
  {

  typename MetricType::Pointer         metric        = MetricType::New();
  typename TransformType::Pointer      transform     = TransformType::New();
  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
  typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  typename FixedImagePyramidType::Pointer fixedImagePyramid =
    FixedImagePyramidType::New();
  typename MovingImagePyramidType::Pointer movingImagePyramid =
    MovingImagePyramidType::New();
  typename RegistrationType::Pointer   registration  = RegistrationType::New();

  /******************************************************************
   * Set up the optimizer.
   ******************************************************************/

  // Instantiate an Observer to report the progress of the Optimization
  typename CommandIterationType::Pointer iterationCommand = CommandIterationType::New();
  iterationCommand->SetOptimizer(  optimizer.GetPointer() );

  // Testing optimizer parameter estimator
  typedef OptimizerParameterEstimator< FixedImageType,
                                        MovingImageType,
                                        TransformType >
                                        OptimizerParameterEstimatorType;
  typename OptimizerParameterEstimatorType::Pointer parameterEstimator
                                        = OptimizerParameterEstimatorType::New();

  parameterEstimator->SetFixedImage(fixedImage);
  parameterEstimator->SetMovingImage(movingImage);
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromShift);
  //parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromJacobian);

  //parameterEstimator->DebugOn();
  //optimizer->DebugOn();

  optimizer->SetOptimizerHelper( parameterEstimator );
  // Setting optimizer parameter done

  // need to minimize for Normalized Correlation
  optimizer->SetMaximize(false);

  /******************************************************************
   * Set up the metric.
   ******************************************************************/
  metric->SetFixedImageRegion( fixedImage->GetBufferedRegion() );

  /******************************************************************
   * Set up the registrator.
   ******************************************************************/

  // connect up the components
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetFixedImage( fixedImage );
  registration->SetMovingImage( movingImage );
  registration->SetInterpolator( interpolator );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );
  registration->SetFixedImageRegion( fixedImage->GetBufferedRegion() );

  // set initial parameters to identity
  typename RegistrationType::ParametersType initialParameters(
    transform->GetNumberOfParameters() );

  initialParameters.Fill( 0.0 );
  //initialize identity transform
  if (dimension == 3)
    {
    initialParameters[0] = 1.0;
    initialParameters[4] = 1.0;
    initialParameters[8] = 1.0;
    }
  else if (dimension == 2)
    {
    initialParameters[0] = 1.0;
    initialParameters[3] = 1.0;
    }
  /******************************************************************
   * Attach registration to a simple UI and run registration
   ******************************************************************/
  SimpleMultiResolutionImageRegistrationUI2<RegistrationType>
    simpleUI( registration );

  //FixedImageType::SizeType size = fixedImage->Get;
  typename FixedImageType::RegionType::SizeType inputSize = fixedImage->GetRequestedRegion().GetSize();
  unsigned int minInputSize = inputSize[0];
  for ( unsigned int dim = 1; dim < FixedImageType::ImageDimension; dim++ )
  {
    if (minInputSize > inputSize[dim])
      {
      minInputSize = inputSize[dim];
      }
  }
  unsigned short numberOfLevels = (int) vcl_ceil( vcl_log(minInputSize/8.0f) / vcl_log(2.0f) ) + 1;
  std::cout << " numberOfLevels = " << numberOfLevels << std::endl;

  itk::Array<unsigned int> niter( numberOfLevels );
  //itk::Array<double>       rates( numberOfLevels );

  niter.Fill(numberOfIterations);

  //now learning rates are estimated
  //rates.Fill(1e-4);

  simpleUI.SetNumberOfIterations( niter );
  //simpleUI.SetLearningRates( rates );

  try
    {
    metric->ReinitializeSeed( 121212 );
    registration->SetNumberOfLevels( numberOfLevels );
    registration->SetInitialTransformParameters( initialParameters );

    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Registration failed" << std::endl;
    std::cout << "Reason " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  /***********************************************************
   * Check the results
   ************************************************************/
  typename RegistrationType::ParametersType solution =
    registration->GetLastTransformParameters();

  std::cout << "Solution is: " << solution << std::endl;

  itk::TransformFileWriter::Pointer writer;
  writer = itk::TransformFileWriter::New();
  writer->SetFileName( argv[4] );
  //transform->SetParameters(registration->GetLastTransformParameters());
  writer->AddTransform( transform );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while saving the transforms" << std::endl;
    std::cerr << excp << std::endl;
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

  }

}
