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
 * Tests are disabled since it requires some image files as input.
 */

//#include "itkDemonsImageToImageObjectMetric.h"
#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkQuasiNewtonObjectOptimizer.h"

#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"
#include "itkDeformationFieldTransform.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegistrationMethodImageSource.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itksys/SystemTools.hxx"

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

using namespace itk;

int itkQuasiNewtonObjectOptimizerTest(int argc, char *argv[])
{

  if( argc >= 2 && (argc < 4 || argc > 7))
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " [fixedImageFile movingImageFile ";
    std::cerr << " outputImageFile] ";
    std::cerr << " [numberOfIterations=10] ";
    std::cerr << " [scalarScale=1] [learningRate=100] " << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << argc << std::endl;
  unsigned int numberOfIterations = 100;
  double scalarScale = 1.0;
  double learningRate = 100;
  if( argc >= 7 )
    {
    numberOfIterations = atoi( argv[4] );
    scalarScale = atof( argv[5] );
    learningRate = atof( argv[6] );
    }

  const unsigned int Dimension = 2;
  typedef double PixelType; //I assume png is unsigned short

  // Fixed Image Type
  typedef itk::Image<PixelType,Dimension>               FixedImageType;

  // Moving Image Type
  typedef itk::Image<PixelType,Dimension>               MovingImageType;

  // Size Type
  typedef MovingImageType::SizeType                 SizeType;


  // ImageSource
  typedef itk::testhelper::ImageRegistrationMethodImageSource<
                                  FixedImageType::PixelType,
                                  MovingImageType::PixelType,
                                  Dimension >         ImageSourceType;

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

  //matching intensity histogram
  typedef HistogramMatchingImageFilter<
                                    MovingImageType,
                                    MovingImageType >   MatchingFilterType;
  MatchingFilterType::Pointer matcher = MatchingFilterType::New();

  matcher->SetInput( movingImage );
  matcher->SetReferenceImage( fixedImage );

  matcher->SetNumberOfHistogramLevels( 256 );
  matcher->SetNumberOfMatchPoints( 10 );
  matcher->ThresholdAtMeanIntensityOn();

  matcher->Update();
  movingImage = matcher->GetOutput();

  //create a deformation field transform
  typedef TranslationTransform<double, Dimension>
                                                  TranslationTransformType;
  TranslationTransformType::Pointer translationTransform =
                                                  TranslationTransformType::New();
  translationTransform->SetIdentity();

  //identity transform for fixed image
  typedef IdentityTransform<double, Dimension>    IdentityTransformType;
  IdentityTransformType::Pointer identityTransform =
                                                  IdentityTransformType::New();
  identityTransform->SetIdentity();

  // The metric
  //  typedef DemonsImageToImageObjectMetric< FixedImageType, MovingImageType >
  //                                                                  MetricType;

  typedef ANTSNeighborhoodCorrelationImageToImageObjectMetric< FixedImageType, MovingImageType>
      MetricType;

  MetricType::Pointer metric = MetricType::New();

  // Assign images and transforms.
  // By not setting a virtual domain image or virtual domain settings,
  // the metric will use the fixed image for the virtual domain.
  metric->SetVirtualDomainImage( const_cast<FixedImageType *>(fixedImage.GetPointer()) );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( identityTransform );
  metric->SetMovingTransform( translationTransform );

  Size<Dimension> radSize;
  radSize.Fill(2);
  metric->SetRadius(radSize);

  //Initialize the metric to prepare for use
  metric->Initialize();

  // Optimizer
  //typedef GradientDescentObjectOptimizer  OptimizerType;
  typedef QuasiNewtonObjectOptimizer  OptimizerType;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetScalarScale( scalarScale );
  optimizer->SetUseScalarScale(true);

  // Instantiate an Observer to report the progress of the Optimization
  typedef itk::CommandIterationUpdate<
                                    OptimizerType >      CommandIterationType;
  CommandIterationType::Pointer iterationCommand = CommandIterationType::New();
  iterationCommand->SetOptimizer(  optimizer.GetPointer() );

  std::cout << "Start optimization..." << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Scalar scale: " << scalarScale << std::endl
            << "Learning rate: " << learningRate << std::endl;
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
    return EXIT_FAILURE;
    }

  std::cout << "...finished. " << std::endl
            << "StopCondition: " << optimizer->GetStopConditionDescription()
            << std::endl
            << "Metric: NumberOfValidPoints: "
            << metric->GetNumberOfValidPoints()
            << std::endl;

  //
  // results
  //
  std::cout << " Results: Transform = " << translationTransform->GetParameters() << std::endl;

  return EXIT_SUCCESS;
}
