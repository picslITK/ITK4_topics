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

#include "itkMattesMutualInformationImageToImageObjectMetric.h"
#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkDemonsImageToImageObjectMetric.h"
#include "itkGradientDescentObjectOptimizer.h"
#include "itkQuasiNewtonObjectOptimizer.h"
#include "itkOptimizerParameterEstimator.h"

#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkEuler2DTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"

#include "itkCastImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itksys/SystemTools.hxx"
#include "itkResampleImageFilter.h"

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

int itkMattesMutualInformationImageToImageObjectRegistrationTest(int argc, char *argv[])
{

  if( argc < 4 || argc > 8)
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile movingImageFile ";
    std::cerr << " outputImageFile ";
    std::cerr << " [numberOfIterations] ";
    std::cerr << " [learningRate] [deformationLearningRate] " << std::endl;
    std::cerr << "For test purpose, return PASSED here." << std::endl;
    std::cout << "Test PASSED." << std::endl;
    return EXIT_SUCCESS;
    }

  std::cout << argc << std::endl;
  unsigned int numberOfIterations = 10;
  double learningRate = 0.1;
  double deformationLearningRate = 1;
  if( argc >= 5 )
    numberOfIterations = atoi( argv[4] );
  if( argc >= 6 )
    learningRate = atof( argv[5] );
  if( argc == 7 )
    deformationLearningRate = atof( argv[6] );
  std::cout << " iterations "<< numberOfIterations << " learningRate "<<learningRate << std::endl;

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

  //get the images
  fixedImageReader->Update();
  FixedImageType::Pointer  fixedImage = fixedImageReader->GetOutput();
  movingImageReader->Update();
  MovingImageType::Pointer movingImage = movingImageReader->GetOutput();

  /** define a resample filter that will ultimately be used to deform the image */
  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();


  /** create a composite transform holder for other transforms  */
  typedef itk::CompositeTransform<double, Dimension>    CompositeType;
  typedef CompositeType::ScalarType                     ScalarType;

  CompositeType::Pointer compositeTransform = CompositeType::New();

  //create an affine transform
  typedef AffineTransform<double, Dimension>
  //  typedef Euler2DTransform<double>
                                                    AffineTransformType;
  AffineTransformType::Pointer affineTransform = AffineTransformType::New();
  affineTransform->SetIdentity();
  std::cout <<" affineTransform params " << affineTransform->GetParameters() << std::endl;
  typedef GaussianSmoothingOnUpdateDisplacementFieldTransform<
                                                    double, Dimension>
                                                     DisplacementTransformType;
  DisplacementTransformType::Pointer displacementTransform =
                                              DisplacementTransformType::New();
  typedef DisplacementTransformType::DisplacementFieldType
                                                         DisplacementFieldType;
  DisplacementFieldType::Pointer field = DisplacementFieldType::New();

  //set the field to be the same as the fixed image region, which will
  // act by default as the virtual domain in this example.
  field->SetRegions( fixedImage->GetLargestPossibleRegion() );
  //make sure the field has the same spatial information as the image
  field->CopyInformation( fixedImage );
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
  displacementTransform->SetGaussianSmoothingSigma( 5 );

  //identity transform for fixed image
  typedef IdentityTransform<double, Dimension> IdentityTransformType;
  IdentityTransformType::Pointer identityTransform =
                                                  IdentityTransformType::New();
  identityTransform->SetIdentity();

  // The metric
  Size<Dimension> radSize;
  radSize.Fill(4);
  typedef MattesMutualInformationImageToImageObjectMetric < FixedImageType, MovingImageType >
  //typedef ANTSNeighborhoodCorrelationImageToImageObjectMetric < FixedImageType, MovingImageType >
  //typedef DemonsImageToImageObjectMetric< FixedImageType, MovingImageType >
                                                                  MetricType;
  MetricType::Pointer metric = MetricType::New();
  metric->SetNumberOfHistogramBins(20);
  //metric->SetRadius(radSize);//antscc

  // Assign images and transforms.
  // By not setting a virtual domain image or virtual domain settings,
  // the metric will use the fixed image for the virtual domain.
  metric->SetVirtualDomainImage( fixedImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( identityTransform );
  compositeTransform->AddTransform( affineTransform );
  compositeTransform->SetAllTransformsToOptimizeOn(); //Set back to optimize all.
  compositeTransform->SetOnlyMostRecentTransformToOptimizeOn(); //set to optimize the displacement field
  metric->SetMovingTransform( compositeTransform );
  bool prewarp = false;
  metric->SetPreWarpMovingImage( prewarp );
  metric->SetPreWarpFixedImage( prewarp );
  bool gaussian = false;
  metric->SetUseMovingGradientRecursiveGaussianImageFilter( gaussian );
  metric->SetUseFixedGradientRecursiveGaussianImageFilter( gaussian );
  metric->Initialize();

  typedef GradientDescentObjectOptimizer  OptimizerType;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->StartOptimization();
  /*
  // Optimizer with parameter estimator
  typedef QuasiNewtonObjectOptimizer  QOptimizerType;
  QOptimizerType::Pointer  qoptimizer = QOptimizerType::New();
  qoptimizer->SetMetric( metric );
  qoptimizer->SetNumberOfIterations( numberOfIterations );
  qoptimizer->SetLineSearchEnabled(false);
  typedef itk::OptimizerParameterEstimator< MetricType > OptimizerParameterEstimatorType;
  OptimizerParameterEstimatorType::Pointer parameterEstimator = OptimizerParameterEstimatorType::New();
  parameterEstimator->SetMetric(metric);
  parameterEstimator->SetTransformForward(true);
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromShift);
  qoptimizer->SetOptimizerParameterEstimator( parameterEstimator );


  std::cout << "Start optimization..." << std::endl
            << "Number of threads: " << metric->GetNumberOfThreads() << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Learning rate: " << learningRate << std::endl
            << "PreWarpImages: " << metric->GetPreWarpMovingImage() << std::endl;
  try
    {
    qoptimizer->StartOptimization();
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
  std::cout << " parameters "  << std::endl;
  std::cout <<  compositeTransform->GetParameters() << std::endl;
*/

  // now add the displacement field to the composite transform
  compositeTransform->AddTransform( displacementTransform );
  compositeTransform->SetOnlyMostRecentTransformToOptimizeOn(); //set to optimize the displacement field
  // Optimizer
  typedef GradientDescentObjectOptimizer  OptimizerType;
  OptimizerType::Pointer  defoptimizer = OptimizerType::New();
  defoptimizer->SetMetric( metric );
  defoptimizer->SetLearningRate( deformationLearningRate );
  defoptimizer->SetNumberOfIterations( numberOfIterations );

  metric->Initialize();

  std::cout << "Start optimization..." << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Deformation learning rate: " << deformationLearningRate << std::endl
            << "PreWarpMovingImage: " << metric->GetPreWarpMovingImage() << std::endl
            << "PreWarpFixedImage: " << metric->GetPreWarpFixedImage() << std::endl
            << "Use_Moving_GradientRecursiveGaussianImageFilter: "
            << metric->GetUseMovingGradientRecursiveGaussianImageFilter()
            << std::endl
            << "Use_Fixed_GradientRecursiveGaussianImageFilter: "
            << metric->GetUseFixedGradientRecursiveGaussianImageFilter()
            << std::endl;
  try
    {
    defoptimizer->StartOptimization();
    }
  catch( ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during deformation Optimization:" << std::endl;
    std::cout << e.GetLocation() << std::endl;
    std::cout << e.GetDescription() << std::endl;
    std::cout << e.what()    << std::endl;
    std::cout << "Test FAILED." << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "...finished. " << std::endl;


  //warp the image with the displacement field
  resample->SetTransform( compositeTransform );
  resample->SetInput( movingImageReader->GetOutput() );
  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 0 );
  resample->Update();
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
  typedef Image< OutputPixelType, Dimension > OutputImageType;
  typedef CastImageFilter<
                        MovingImageType,
                        OutputImageType >     CastFilterType;
  typedef ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();
  writer->SetFileName( argv[3] );
  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput() );
  writer->Update();

  std::cout << "Test PASSED." << affineTransform->GetParameters() << std::endl;
  return EXIT_SUCCESS;

}
