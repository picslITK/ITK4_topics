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
 * Test program for ANTSNeighborhoodCorrelationImageToImageObjectMetric and
 * GradientDescentObjectOptimizer classes, using a pair of input images.
 *
 */

//#include "itkDemonsImageToImageObjectMetric.h"
#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkGradientDescentObjectOptimizer.h"

#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itksys/SystemTools.hxx"

//We need this as long as we have to define ImageToData as a fwd-declare
// in itkImageToImageObjectMetric.h
#include "itkImageToData.h"

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

int itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest(int argc, char *argv[])
{

  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile movingImageFile ";
    std::cerr << " outputImageFile ";
    std::cerr << " [numberOfIterations=100] ";
    std::cerr << " [learningRate=100] " << std::endl;
    std::cerr << " [usePreWarp=1 | 0]" << std::endl;
    std::cerr << "For test purpose, return PASSED here." << std::endl;
    std::cout << "Test PASSED." << std::endl;
    return EXIT_SUCCESS;
    }
//  std::cout << argc << std::endl;
  unsigned int numberOfIterations = 100;
  double learningRate = 100;
  bool preWarp = true;
  if( argc >= 5 )
    numberOfIterations = atoi( argv[4] );
  if( argc >= 6 )
    learningRate = atof( argv[5] );
  if ( argc >= 7 )
    preWarp = (atoi(argv[6]) != 0);

  const unsigned int Dimension = 2;
  typedef double PixelType; //I assume png is unsigned short

  typedef Image< PixelType, Dimension >  FixedImageType;
  typedef Image< PixelType, Dimension >  MovingImageType;

  typedef ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer fixedImageReader   = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName( argv[1] );
  movingImageReader->SetFileName( argv[2] );

  //matching intensity histogram
  typedef HistogramMatchingImageFilter<
                                    MovingImageType,
                                    MovingImageType >   MatchingFilterType;
  MatchingFilterType::Pointer matcher = MatchingFilterType::New();

  matcher->SetInput( movingImageReader->GetOutput() );
  matcher->SetReferenceImage( fixedImageReader->GetOutput() );

  matcher->SetNumberOfHistogramLevels( 256 );
  matcher->SetNumberOfMatchPoints( 10 );
  matcher->ThresholdAtMeanIntensityOn();
  //get the images
  fixedImageReader->Update();
  FixedImageType::Pointer  fixedImage = fixedImageReader->GetOutput();
  movingImageReader->Update();
  matcher->Update();
  MovingImageType::Pointer movingImage = matcher->GetOutput();
  // MovingImageType::Pointer movingImage = movingImageReader->GetOutput();


  //create a displacement field transform
  typedef TranslationTransform<double, Dimension>
                                                    TranslationTransformType;
  TranslationTransformType::Pointer translationTransform =
                                                  TranslationTransformType::New();
//  translationTransform->SetIdentity();
  TranslationTransformType::ParametersType initial_translation;
  initial_translation.SetSize(translationTransform->GetNumberOfParameters());
  initial_translation[0] = 3;
  initial_translation[1] = -5;

  translationTransform->SetParameters(initial_translation);

  typedef GaussianSmoothingOnUpdateDisplacementFieldTransform<double, Dimension>
                                                    DisplacementTransformType;
  DisplacementTransformType::Pointer displacementTransform =
                                              DisplacementTransformType::New();
  typedef DisplacementTransformType::DisplacementFieldType DisplacementFieldType;
  DisplacementFieldType::Pointer field = DisplacementFieldType::New();

  //set the field to be the same as the fixed image region, which will
  // act by default as the virtual domain in this example.
  field->SetRegions( fixedImage->GetLargestPossibleRegion() );
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


  typedef ANTSNeighborhoodCorrelationImageToImageObjectMetric< FixedImageType, MovingImageType>
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
  // metric->SetMovingTransform( translationTransform );

  Size<Dimension> radSize;
  radSize.Fill(2);
  metric->SetRadius(radSize);


  metric->SetPreWarpMovingImage( preWarp );
  metric->SetUseMovingGradientRecursiveGaussianImageFilter( false );


  //Initialize the metric to prepare for use
  metric->Initialize();

  // Optimizer
  typedef GradientDescentObjectOptimizer  OptimizerType;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetNumberOfIterations( numberOfIterations );

  std::cout << "Start optimization..." << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Learning rate: " << learningRate << std::endl
            << "CC radius: " << metric->GetRadius() << std::endl
            << "CC prewarp: " << metric->GetPreWarpMovingImage() << std::endl
            << "CC number of threads: " << metric->GetNumberOfThreads() << std::endl;


//  std::cout << "initial para: " << translationTransform->GetParameters() << std::endl;


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

//  std::cout << "final para: " << translationTransform->GetParameters() << std::endl;

  field = displacementTransform->GetDisplacementField();
  std::cout << "LargestPossibleRegion: " << field->GetLargestPossibleRegion()
            << std::endl;
  ImageRegionIteratorWithIndex< DisplacementFieldType > it( field, field->GetLargestPossibleRegion() );

  //warp the image with the displacement field
  typedef WarpImageFilter<
                          MovingImageType,
                          MovingImageType,
                          DisplacementFieldType  >     WarperType;
  typedef LinearInterpolateImageFunction<
                                   MovingImageType,
                                   double          >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  WarperType::Pointer warper = WarperType::New();
  warper->SetInput( movingImage );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixedImage->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOrigin() );
  warper->SetOutputDirection( fixedImage->GetDirection() );

  warper->SetDisplacementField( displacementTransform->GetDisplacementField() );

  //write out the displacement field
  typedef ImageFileWriter< DisplacementFieldType >  DisplacementWriterType;
  DisplacementWriterType::Pointer      displacementwriter =  DisplacementWriterType::New();
  std::string outfilename( argv[3] );
  std::string ext = itksys::SystemTools::GetFilenameExtension( outfilename );
  std::string defout=outfilename + std::string("_def") + ext;
  displacementwriter->SetFileName( defout.c_str() );
  displacementwriter->SetInput( displacementTransform->GetDisplacementField() );
  displacementwriter->Update();

  //write the warped image into a file
  typedef double                              OutputPixelType;
  typedef Image< OutputPixelType, Dimension > OutputImageType;
  typedef CastImageFilter<
                        MovingImageType,
                        OutputImageType >     CastFilterType;
  typedef ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  caster->SetInput( warper->GetOutput() );
  writer->SetInput( caster->GetOutput() );
  writer->Update();

  std::cout << "Test PASSED." << std::endl;
  return EXIT_SUCCESS;
}
