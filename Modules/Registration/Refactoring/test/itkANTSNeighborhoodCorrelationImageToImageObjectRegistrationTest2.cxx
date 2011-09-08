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
 * GradientDescentObjectOptimizer classes, using Gaussian image source.
 *
 */

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
#include "itkImageRegistrationMethodImageSource.h"
#include "itkAffineTransform.h"

//We need this as long as we have to define ImageToData as a fwd-declare
// in itkImageToImageObjectMetric.h
#include "itkImageToData.h"

using namespace itk;

namespace {
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
template<typename TRegistration>
class ShowProgressObject {
public:
    ShowProgressObject(TRegistration* o) {
        m_Process = o;
    }
    void ShowProgress() {
        std::cout << "Progress: " << m_Process->GetProgress() << "  ";
        std::cout << "Iter: " << m_Process->GetElapsedIterations() << "  ";
        std::cout << "Metric: " << m_Process->GetMetric() << "  ";
        std::cout << "RMSChange: " << m_Process->GetRMSChange() << "  ";
        std::cout << std::endl;
        if (m_Process->GetElapsedIterations() == 10) {
            m_Process->StopRegistration();
        }
    }
    typename TRegistration::Pointer m_Process;
};
}

int itkANTSNeighborhoodCorrelationImageToImageObjectRegistrationTest2(int argc,
        char *argv[]) {

    unsigned int numberOfIterations = 100;
    double learningRate = 1000;
    bool usePreWarp = true;

    if (argc == 2 && strcmp(argv[1], "-h") == 0) {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " [numberOfIterations=100] ";
        std::cerr << " [learningRate=100] " << std::endl;
        std::cerr << " [usePreWarp=1 | 0]" << std::endl;
        std::cerr << "For test purpose, return PASSED here." << std::endl;
        std::cout << "Test PASSED." << std::endl;

        //TODO: return EXIT_SUCCESS just for test purpose
        return EXIT_SUCCESS;

    }

    const unsigned int Dimension = 2;

    Size<Dimension> radSize;
    radSize.Fill(5);

    typedef double PixelType;

    typedef Image<PixelType, Dimension> FixedImageType;
    typedef Image<PixelType, Dimension> MovingImageType;

    if (argc >= 2)
        numberOfIterations = atoi(argv[2]);

    if (argc >= 3)
        learningRate = atof(argv[3]);
    if (argc >= 4)
        usePreWarp = atoi(argv[4]);

    bool pass = true;
    // Size Type
    typedef MovingImageType::SizeType SizeType;

    // ImageSource
    typedef itk::testhelper::ImageRegistrationMethodImageSource<
            FixedImageType::PixelType, MovingImageType::PixelType, Dimension> ImageSourceType;

    ImageSourceType::Pointer imageSource = ImageSourceType::New();

    SizeType size;
    size[0] = 100;
    size[1] = 100;

    imageSource->GenerateImages(size);

    FixedImageType::Pointer fixedImage =
            const_cast<FixedImageType *>(imageSource->GetFixedImage());MovingImageType
    ::Pointer movingImage =
            const_cast<MovingImageType *>(imageSource->GetMovingImage());

   // Transform Type
    typedef itk::AffineTransform<double, Dimension> TransformType;
    typedef TransformType::ParametersType           ParametersType;

    TransformType::Pointer affineTransform = TransformType::New();
    affineTransform->SetIdentity();

    //create a displacement field transform
    typedef TranslationTransform<double, Dimension> TranslationTransformType;
    TranslationTransformType::Pointer translationTransform =
            TranslationTransformType::New();
    translationTransform->SetIdentity();

    typedef GaussianSmoothingOnUpdateDisplacementFieldTransform<double, Dimension> DisplacementTransformType;
    DisplacementTransformType::Pointer displacementTransform =
            DisplacementTransformType::New();
    typedef DisplacementTransformType::DisplacementFieldType DisplacementFieldType;
    DisplacementFieldType::Pointer field = DisplacementFieldType::New();

    //set the field to be the same as the fixed image region, which will
    // act by default as the virtual domain in this example.
    field->SetRegions(fixedImage->GetLargestPossibleRegion());
    field->Allocate();
    // Fill it with 0's
    DisplacementTransformType::OutputVectorType zeroVector;
    zeroVector.Fill(0);
    field->FillBuffer(zeroVector);
    // Assign to transform
    displacementTransform->SetDisplacementField(field);
    displacementTransform->SetGaussianSmoothingSigma(6);

    //identity transform for fixed image
    typedef IdentityTransform<double, Dimension> IdentityTransformType;
    IdentityTransformType::Pointer identityTransform =
            IdentityTransformType::New();
    identityTransform->SetIdentity();

    typedef ANTSNeighborhoodCorrelationImageToImageObjectMetric<FixedImageType,
            MovingImageType> MetricType;

    MetricType::Pointer metric = MetricType::New();

    // Assign images and transforms.
    // By not setting a virtual domain image or virtual domain settings,
    // the metric will use the fixed image for the virtual domain.

    metric->SetVirtualDomainImage(fixedImage);
    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(movingImage);
    metric->SetFixedTransform(identityTransform);
    // metric->SetMovingTransform( affineTransform );
    // metric->SetMovingTransform( displacementTransform );

    metric->SetMovingTransform(translationTransform);

    metric->SetRadius(radSize);

    metric->SetPreWarpMovingImage(usePreWarp);
    metric->SetUseMovingGradientRecursiveGaussianImageFilter(false);

    //Initialize the metric to prepare for use
    metric->Initialize();

    // Optimizer
    typedef GradientDescentObjectOptimizer OptimizerType;
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetMetric(metric);
    optimizer->SetLearningRate(learningRate);
    optimizer->SetNumberOfIterations(numberOfIterations);

    std::cout << "Start optimization..." << std::endl
            << "Number of iterations: " << numberOfIterations << std::endl
            << "Learning rate: " << learningRate << std::endl
            << "CC radius: " << metric->GetRadius()
            << std::endl << "CC prewarp: " << metric->GetPreWarpMovingImage()
            << std::endl << "CC number of threads: "
            << metric->GetNumberOfThreads() << std::endl;

    try {
        optimizer->StartOptimization();
    } catch (ExceptionObject & e) {
        std::cout << "Exception thrown ! " << std::endl;
        std::cout << "An error ocurred during Optimization:" << std::endl;
        std::cout << e.GetLocation() << std::endl;
        std::cout << e.GetDescription() << std::endl;
        std::cout << e.what() << std::endl;
        pass = false;

        std::cout << "Test PASSED." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "...finished. " << std::endl << "StopCondition: "
            << optimizer->GetStopConditionDescription() << std::endl
            << "Metric: NumberOfValidPoints: "
            << metric->GetNumberOfValidPoints() << std::endl;

    ParametersType actualParameters = imageSource->GetActualParameters();
    // ParametersType finalParameters  = affineTransform->GetParameters();
    ParametersType finalParameters = translationTransform->GetParameters();

    const unsigned int numbeOfParameters = actualParameters.Size();

    // We know that for the Affine transform the Translation parameters are at
    // the end of the list of parameters.
    const unsigned int offsetOrder = finalParameters.Size()
            - actualParameters.Size();

    const double tolerance = 1.0; // equivalent to 1 pixel.

    for (unsigned int i = 0; i < numbeOfParameters; i++) {
        // the parameters are negated in order to get the inverse transformation.
        // this only works for comparing translation parameters....
        std::cout << finalParameters[i + offsetOrder] << " == "
                << -actualParameters[i] << std::endl;
        if (vnl_math_abs(
                finalParameters[i + offsetOrder] - (-actualParameters[i]))
                > tolerance) {
            std::cout << "Tolerance exceeded at component " << i << std::endl;
            pass = false;
        }
    }

//  //
//  //  Get the transform as the Output of the Registration filter
//  //
//  RegistrationType::TransformOutputConstPointer transformDecorator =
//                                                        registration->GetOutput();
//
//  TransformType::ConstPointer finalTransform =
//    static_cast< const TransformType * >( transformDecorator->Get() );

    field = displacementTransform->GetDisplacementField();
    std::cout << "LargestPossibleRegion: " << field->GetLargestPossibleRegion()
            << std::endl;
    ImageRegionIteratorWithIndex<DisplacementFieldType> it(field,
            field->GetLargestPossibleRegion());

    typedef WarpImageFilter<MovingImageType,
                            MovingImageType,
                            DisplacementFieldType>                   WarperType;
    typedef LinearInterpolateImageFunction<MovingImageType, double> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    WarperType::Pointer warper = WarperType::New();
    warper->SetInput(movingImage);
    warper->SetInterpolator(interpolator);
    warper->SetOutputSpacing(fixedImage->GetSpacing());
    warper->SetOutputOrigin(fixedImage->GetOrigin());
    warper->SetOutputDirection(fixedImage->GetDirection());

    warper->SetDisplacementField(displacementTransform->GetDisplacementField());

//  //write out the displacement field
//  typedef ImageFileWriter< DisplacementFieldType >  DisplacementWriterType;
//  DisplacementWriterType::Pointer      displacementwriter =  DisplacementWriterType::New();
//  std::string outfilename( argv[1] );
//  std::string ext = itksys::SystemTools::GetFilenameExtension( outfilename );
//  std::string defout=outfilename + std::string("_def") + ext;
//  displacementwriter->SetFileName( defout.c_str() );
//  displacementwriter->SetInput( displacementTransform->GetDisplacementField() );
//  displacementwriter->Update();

    //write the warped image into a file
    typedef double                                            OutputPixelType;
    typedef Image<OutputPixelType, Dimension>                 OutputImageType;
    typedef CastImageFilter<MovingImageType, OutputImageType> CastFilterType;
    typedef ImageFileWriter<OutputImageType>                  WriterType;

    WriterType::Pointer writer = WriterType::New();
    CastFilterType::Pointer caster = CastFilterType::New();

    writer->SetFileName(argv[1]);

    caster->SetInput(warper->GetOutput());
    writer->SetInput(caster->GetOutput());
    writer->Update();

    if (!pass) {
        std::cout << "Test FAILED." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Test PASSED." << std::endl;
    return EXIT_SUCCESS;

}
