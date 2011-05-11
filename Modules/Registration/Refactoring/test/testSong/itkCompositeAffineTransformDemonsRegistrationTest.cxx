/*
 * itkCompositeAffineTransformDemonsRegistrationTest.cxx
 *
 * Example: composite transform and multi-threaded registration
 *
 *  transform: centered composite affine
 *  metric: demons metric
 *  optmizer: multithreaded optimizer
 *
 *  Created on: May 11, 2011
 *      Author: songgang
 */

#include "itkImage.h"
#include "itkCenteredCompositeTransform.h"
#include "itkDemonsImageToImageMetric.h"
#include "itkObjectToObjectThreadedMetricOptimizer.h"

#include "itkImageFileReader.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "itkCompositeTransform.h"
#include "itkCenteredCompositeTransform.h"
#include <itkTimeProbe.h>
#include "vcl_cmath.h"
#include "itkScaleTransform.h"
#include "itkTranslationTransform.h"
#include "itkShear2DTransform.h"
#include "itkRotate2DTransform.h"

#include "itkIdentityTransform.h"

#include "itkObjectToObjectThreadedMetricOptimizer.h"
#include "itkImageToData.h"

template<class ImagePointerType>
void CompositeAffineTransformDemonsRegistrationTest(
        const ImagePointerType & fixImage,
        const ImagePointerType & movImage);



int main(int argc, char** argv){

    const int Dimension = 2;
    typedef itk::Image<float, Dimension> ImageType;
    typedef ImageType::Pointer ImagePointerType;

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    ImageReaderType::Pointer imageReader = ImageReaderType::New();


    imageReader->SetFileName(argv[0]);
    imageReader->Update();
    ImagePointerType fixImage = imageReader->GetOutput();

    imageReader->SetFileName(argv[1]);
    imageReader->Update();
    ImagePointerType movImage = imageReader->GetOutput();

    CompositeAffineTransformDemonsRegistrationTest(fixImage, movImage);



    return EXIT_SUCCESS;
}


template<class CompositeTransformPointerType>
void CreateCompositeAffineTransform(CompositeTransformPointerType &composite_transform) {

    typedef typename CompositeTransformPointerType::ObjectType CompositeTransformType;
    const int Dim = CompositeTransformType::InputDimension;


    typedef itk::ScaleTransform<double, Dim> ScaleTransformType;
    typedef itk::Shear2DTransform<double, Dim> Shear2DTransformType;
    typedef itk::Rotate2DTransform<double> Rotation2DTransformType;
    typedef itk::TranslationTransform<double, Dim> TranslationTransformType;

    typename ScaleTransformType::Pointer scale_transform = ScaleTransformType::New();
    scale_transform->SetIdentity();
    typename ScaleTransformType::ParametersType s1(Dim);
    for(int i=0; i<Dim; i++) s1[i]=1; //(rand() % 100) / 50.0 ; //
    scale_transform->SetParameters(s1);


    typename Shear2DTransformType::Pointer shear_transform = Shear2DTransformType::New();
    typename Shear2DTransformType::ParametersType k1(1);
    for(int i=0; i<1; i++) k1[i]=0; //(rand() % 100) / 50.0 - 1.0 ; //
    shear_transform->SetParameters(k1);


    typename Rotation2DTransformType::Pointer rotation_transform = Rotation2DTransformType::New();
    typename Rotation2DTransformType::ParametersType r1(1);
    r1.Fill(0); // measured in pi
    rotation_transform->SetParameters(r1);


    typename TranslationTransformType::Pointer translation_transform = TranslationTransformType::New();
    typename TranslationTransformType::ParametersType t1(Dim);
    for(int i=0; i<Dim; i++) t1[i]=0; // (rand() % 100) / 50.0 - 1.0 ; //
    translation_transform->SetParameters(t1);

    std::cout << "r.center=" << rotation_transform->GetCenter() << std::endl;
    std::cout << "r.fixed=" << rotation_transform->GetFixedParameters() << std::endl;
    std::cout << "r.trans=" << rotation_transform->GetTranslation() << std::endl;


    typename CompositeTransformType::Pointer comp = CompositeTransformType::New();

    // add in the reverse order of applying to the point
    comp->AddTransform(translation_transform);
    comp->AddTransform(rotation_transform);
    comp->AddTransform(scale_transform);
    comp->AddTransform(shear_transform);
}



template<class ImagePointerType>
void CompositeAffineTransformDemonsRegistrationTest(
        const ImagePointerType & fixImage,
        const ImagePointerType & movImage) {

    typedef typename ImagePointerType::ObjectType ImageType;
    const int Dim = ImageType::ImageDimension;

    std::cout << "fixImage: " << fixImage << std::endl;
    std::cout << "movImage: " << movImage << std::endl;

    typedef itk::CenteredCompositeTransform<double, Dim> CompositeTransformType;
    typename CompositeTransformType::Pointer comp = CompositeTransformType::New();
    CreateCompositeAffineTransform(comp);


    typedef itk::IdentityTransform<double, Dim> IdentityTransformType;
    typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();


    typedef itk::DemonsImageToImageMetric<ImageType, ImageType> MetricType;
    typename MetricType::Pointer metric = MetricType::New();

    metric->SetFixedImage(fixImage);
    metric->SetMovingImage(movImage);
    metric->SetFixedImageTransform(identityTransform);
    metric->SetMovingImageTransform(comp);
    metric->SetVirtualDomainSize(fixImage->GetRequestedRegion().GetSize());
    metric->SetVirtualDomainIndex(fixImage->GetRequestedRegion().GetIndex());
    metric->SetVirtualDomainSpacing(fixImage->GetSpacing());
    metric->SetVirtualDomainOrigin(fixImage->GetOrigin());
    metric->SetVirtualDomainDirection(fixImage->GetDirection());
    // Compute one iteration of the metric
    metric->Initialize();



    int numberOfThreads = 4;

    typedef itk::ObjectToObjectThreadedMetricOptimizer<MetricType> OptimizerType;
    typedef itk::ImageToData<Dim, OptimizerType> ThreaderType;

    // pseudo code
    typename ThreaderType::Pointer metricThreader = ThreaderType::New();
    OptimizerType metricOptimizer;

    metricOptimizer.metric = metric;

    typename ImageType::RegionType inboundary_region = fixImage->GetRequestedRegion();

    metricThreader->SetNumberOfThreads(numberOfThreads);
    metricThreader->m_OverallRegion = inboundary_region ;
    metricThreader->m_Holder = &metricOptimizer;
    metricThreader->ThreadedGenerateData = OptimizerType::ComputeMetricValueInRegionThreaded;
    metricOptimizer.BeforeThreadedGenerateData(numberOfThreads);
    metricThreader->GenerateData();


    float energy = static_cast<float> (metricOptimizer.AccumulateMeasuresFromAllThreads());

    std::cout << "metric = " << energy << std::endl;
    return;

}
