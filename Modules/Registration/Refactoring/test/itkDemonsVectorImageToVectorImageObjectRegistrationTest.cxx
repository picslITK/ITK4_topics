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

/**
 * Test program for DemonImageToImageObjectMetric and
 * GradientDescentObjectOptimizer classes.
 *
 * Perform a registration using user-supplied images.
 * No numerical verification is performed. Test passes as long
 * as no exception occurs.
 */

#include "itkDemonsVectorImageToVectorImageObjectMetric.h"
#include "itkGradientDescentObjectOptimizer.h"

#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"

#include "itkHistogramMatchingImageFilter.h"
//#include "itkCastImageFilter.h"
//#include "itkWarpImageFilter.h"
//These two are needed as long as we're using fwd-declarations in
//DisplacementFieldTransfor:
#include "itkVectorInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itksys/SystemTools.hxx"

// #include "itkMinimumMaximumImageCalculator.h"

using namespace itk;

namespace{
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
template<typename TRegistration>
class ShowProgressObject
{
public:
  ShowProgressObject(TRegistration* o)
    {m_Process = o;}
  void ShowProgress()
    {
    std::cout << "Progress: " << m_Process->GetProgress() << "  ";
    std::cout << "Iter: " << m_Process->GetElapsedIterations() << "  ";
    std::cout << "Metric: "   << m_Process->GetMetric()   << "  ";
    std::cout << "RMSChange: " << m_Process->GetRMSChange() << "  ";
    std::cout << std::endl;
    if ( m_Process->GetElapsedIterations() == 10 )
      { m_Process->StopRegistration(); }
    }
  typename TRegistration::Pointer m_Process;
};
}

int itkDemonsVectorImageToVectorImageObjectRegistrationTest(int argc, char *argv[])
{

  if( argc < 4 || argc > 7)
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile movingImageFile ";
    std::cerr << " outputImageFile ";
    std::cerr << " [numberOfIterations] ";
    std::cerr << " [learningRate] " << std::endl;
    std::cerr << "For test purpose, return PASSED here." << std::endl;
    std::cout << "Test PASSED." << std::endl;
    return EXIT_SUCCESS;
    }

  std::cout << argc << std::endl;
  unsigned int numberOfIterations = 10;
  double learningRate = 0.1;
  if( argc >= 5 )
    numberOfIterations = atoi( argv[4] );
  if( argc == 6 )
    learningRate = atof( argv[5] );

  const unsigned int Dimension = 2;
  typedef double PixelType; //I assume png is unsigned short

  typedef VectorImage< PixelType, Dimension >  FixedImageType;
  typedef VectorImage< PixelType, Dimension >  MovingImageType;

  typedef ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer fixedImageReader   = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName( argv[1] );
  movingImageReader->SetFileName( argv[2] );

  //matching intensity histogram
  //typedef HistogramMatchingImageFilter<
  //                                  MovingImageType,
  //                                  MovingImageType >   MatchingFilterType;
  //MatchingFilterType::Pointer matcher = MatchingFilterType::New();

  //matcher->SetInput( movingImageReader->GetOutput() );
  //matcher->SetReferenceImage( fixedImageReader->GetOutput() );

  //matcher->SetNumberOfHistogramLevels( 256 );
  //matcher->SetNumberOfMatchPoints( 10 );
  //matcher->ThresholdAtMeanIntensityOn();
  //get the images
  fixedImageReader->Update();
  FixedImageType::Pointer  fixedImage = fixedImageReader->GetOutput();
  movingImageReader->Update();
  //matcher->Update();
  //MovingImageType::Pointer movingImage = matcher->GetOutput();
  MovingImageType::Pointer movingImage = movingImageReader->GetOutput();

  //create a displacement field transform
  typedef TranslationTransform<double, Dimension>
                                                    TranslationTransformType;
  TranslationTransformType::Pointer translationTransform =
                                                  TranslationTransformType::New();
  translationTransform->SetIdentity();

  typedef GaussianSmoothingOnUpdateDisplacementFieldTransform<double,
                                                              Dimension>
                                                     DisplacementTransformType;
  DisplacementTransformType::Pointer displacementTransform =
                                              DisplacementTransformType::New();
  typedef DisplacementTransformType::DisplacementFieldType
                                                         DisplacementFieldType;
  DisplacementFieldType::Pointer field = DisplacementFieldType::New();

  //set the field to be the same as the fixed image region, which will
  // act by default as the virtual domain in this example.
  field->SetRegions( fixedImage->GetLargestPossibleRegion() );
  std::cout << "fixedImage->GetLargestPossibleRegion(): "
            << fixedImage->GetLargestPossibleRegion() << std::endl
            << "fixedImage->GetBufferedRegion(): "
            << fixedImage->GetBufferedRegion() << std::endl;
  field->Allocate();
  // Fill it with 0's
  DisplacementTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );
  // Assign to transform
  displacementTransform->SetDisplacementField( field );
  displacementTransform->SetGaussianSmoothingSigma( 6 );

  //identity transform for fixed image
  typedef IdentityTransform<double, Dimension> IdentityTransformType;
  IdentityTransformType::Pointer identityTransform =
                                                  IdentityTransformType::New();
  identityTransform->SetIdentity();

  // The metric
  typedef DemonsVectorImageToVectorImageObjectMetric< FixedImageType, MovingImageType >
                                                                  MetricType;
  MetricType::Pointer metric = MetricType::New();

  // Assign images and transforms.
  // By not setting a virtual domain image or virtual domain settings,
  // the metric will use the fixed image for the virtual domain.
  metric->SetVirtualDomainImage( fixedImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( identityTransform );
  metric->SetMovingTransform( displacementTransform );
  //  metric->SetMovingTransform( translationTransform );

  metric->SetPreWarpImages( true );
  metric->SetPrecomputeImageGradient( ! metric->GetPreWarpImages() );

  //Initialize the metric to prepare for use
  metric->Initialize();

  MetricType::MeasureType valueReturn;
  MetricType::DerivativeType derivativeReturn;
  metric->GetValueAndDerivative( valueReturn, derivativeReturn );
  std::cout << "Initial metric: " << valueReturn << std::endl;

  // Optimizer
  typedef GradientDescentObjectOptimizer  OptimizerType;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetNumberOfIterations( numberOfIterations );

  std::cout << "Start optimization..." << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Learning rate: " << learningRate << std::endl
            << "PreWarpImages: " << metric->GetPreWarpImages() << std::endl;
  try
    {
    optimizer->StartOptimization();
    }
  catch( ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization:" << std::endl;
    std::cout << e.GetLocation() << std::endl;
    std::cout << e.GetDescription() << std::endl;
    std::cout << e.what()    << std::endl;
    std::cout << "Test FAILED." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "...finished. " << std::endl
            << "StopCondition: " << optimizer->GetStopConditionDescription()
            << std::endl
            << "Metric: NumberOfValidPoints: "
            << metric->GetNumberOfValidPoints()
            << std::endl;

  field = displacementTransform->GetDisplacementField();
  std::cout << "LargestPossibleRegion: " << field->GetLargestPossibleRegion()
            << std::endl;
  ImageRegionIteratorWithIndex< DisplacementFieldType > it( field, field->GetLargestPossibleRegion() );
  /* print out a few displacement field vectors */
  /*std::cout
      << "First few elements of first few rows of final displacement field:"
      << std::endl;
  for(unsigned int i=0; i< 5; i++ )
    {
     for(unsigned int j=0; j< 5; j++ )
      {
      DisplacementFieldType::IndexType index;
      index[0] = i;
      index[1] = j;
      it.SetIndex(index);
      std::cout << it.Value() << " ";
      }
    std::cout << std::endl;
    }
  */

  //
  // results
  //
  //  std::cout << " result " << translationTransform->GetParameters() << std::endl;
  //warp the image with the displacement field
  typedef VectorResampleImageFilter<
                          MovingImageType,
                          MovingImageType,
                          double  >     WarperType;
  typedef VectorLinearInterpolateImageFunction<
                                   MovingImageType,
                                   double          >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  WarperType::Pointer warper = WarperType::New();
  warper->SetInput( movingImage );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixedImage->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOrigin() );
  warper->SetOutputDirection( fixedImage->GetDirection() );
  warper->SetTransform( displacementTransform );
  warper->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  //warper->SetDisplacementField( displacementTransform->GetDisplacementField() );

  //write out the displacement field
  typedef ImageFileWriter< DisplacementFieldType >  DisplacementWriterType;
  DisplacementWriterType::Pointer      displacementwriter =  DisplacementWriterType::New();
  std::string outfilename( argv[3] );
  std::string ext = itksys::SystemTools::GetFilenameExtension( outfilename );
  std::string name = itksys::SystemTools::GetFilenameWithoutExtension( outfilename );
  std::string defout = name + std::string("_def") + ext;
  displacementwriter->SetFileName( defout.c_str() );
  displacementwriter->SetInput( displacementTransform->GetDisplacementField() );
  displacementwriter->Update();

  //write the warped image into a file
  typedef double                              OutputPixelType;
  //typedef Image< OutputPixelType, Dimension > OutputImageType;
  //typedef CastImageFilter<
  //                      MovingImageType,
  //                      OutputImageType >     CastFilterType;
  typedef ImageFileWriter< MovingImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  //CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  //caster->SetInput( warper->GetOutput() );
  writer->SetInput( warper->GetOutput() );

  writer->Update();

  std::cout << "Test PASSED." << std::endl;
  return EXIT_SUCCESS;

}
