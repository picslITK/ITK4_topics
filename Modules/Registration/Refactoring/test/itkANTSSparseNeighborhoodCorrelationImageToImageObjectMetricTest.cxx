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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkDemonsImageToImageObjectMetric.h"
#include "itkIdentityTransform.h"
#include "itkDeformationFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"

#include "itkImageRegionConstIterator.h"

#include "itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric.h"

using namespace itk;

template<class ImagePointerType, class DerivativeType>
void PrintDerivativeAsVectorImage(ImagePointerType image, DerivativeType &derivative, unsigned int vecdim){

    // typedef typename ImageType::ConstPointer ImageConstPointerType;
    // ImageConstPointerType image = imageP;

    typedef typename ImagePointerType::ObjectType ImageType;
    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    unsigned int dim0 = imageRegion.GetSize()[0];
    unsigned int dim1 = imageRegion.GetSize()[1];

    typedef ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.Begin();
    unsigned int cnt = 0;
    for (unsigned int ycnt = 0; ycnt < dim1; ycnt++) {
        for (unsigned int xcnt = 0; xcnt < dim0; xcnt++) {

            std::cout << '[';
            for(unsigned int d = 0; d < vecdim-1; d++){
                std::cout << derivative[cnt*vecdim + d] <<",";
            }
            std::cout << derivative[cnt * vecdim + vecdim-1] << ']' << "\t";

            ++it;
            ++cnt;

        }
        std::cout << std::endl;
    }

    return;

}


template<class ImageType>
void PrintImage(ImageType *imageP) {

    typedef typename ImageType::ConstPointer ImageConstPointerType;
    ImageConstPointerType image = imageP;

    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    unsigned int dim0 = imageRegion.GetSize()[0];
    unsigned int dim1 = imageRegion.GetSize()[1];

    typedef ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.Begin();
//    unsigned int cnt = 0;
    for (unsigned ycnt = 0; ycnt < dim1; ycnt++) {
        for (unsigned xcnt = 0; xcnt < dim0; xcnt++) {
            std::cout << it.Get() << "\t";
            ++it;
        }
        std::cout << std::endl;
    }

    return;
}

template<class ImagePointerType>
void PrintImage(const ImagePointerType &image) {

    typedef typename ImagePointerType::ObjectType ImageType;
    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    unsigned int dim0 = imageRegion.GetSize()[0];
    unsigned int dim1 = imageRegion.GetSize()[1];

    typedef ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.Begin();
//    unsigned int cnt = 0;
    for (unsigned ycnt = 0; ycnt < dim1; ycnt++) {
        for (unsigned xcnt = 0; xcnt < dim0; xcnt++) {
            std::cout << it.Get() << "\t";
            ++it;
        }
        std::cout << std::endl;
    }

    return;
}

int itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetricTest(
        int argc, char * argv[]) {


    MultiThreader::SetGlobalMaximumNumberOfThreads(1);

//  if (argc <= 3) {
//    std::cout << "Args: 2D_fixed_image 2D_moving_image number_of_threads" << std::endl;
//    exit(-1);
//  }

//  const char * filename1 = argv[1];
//  const char * filename2 = argv[2];
//  int number_of_threads = atoi(argv[3]);

    const int ImageDimension = 2;

    typedef itk::Image<double, ImageDimension> ImageType;
    typedef ImageType::Pointer ImagePointerType;
    typedef ImageType::RegionType RegionType;

    typedef itk::Vector<double, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension> VectorImageType;

    typedef itk::Transform<double, ImageDimension> TransformType;
    typedef itk::IdentityTransform<double, ImageDimension> IdentityTransformType;
    typedef itk::CompositeTransform<double, ImageDimension> CompositeTransformType;
    typedef itk::TranslationTransform<double, ImageDimension> TranslationTransformType;
    typedef itk::DeformationFieldTransform<double, ImageDimension> DeformationTransformType;
    typedef DeformationTransformType::DeformationFieldType FieldType;

    IdentityTransformType::Pointer transformFId = IdentityTransformType::New();

    IdentityTransformType::Pointer transformMId = IdentityTransformType::New();
    DeformationTransformType::Pointer transformMdeformation =
            DeformationTransformType::New();
    TranslationTransformType::Pointer transformMtranslation =
            TranslationTransformType::New();
    TranslationTransformType::Pointer transformMtranslation2 =
            TranslationTransformType::New();
    CompositeTransformType::Pointer transformMComp =
            CompositeTransformType::New();
    CompositeTransformType::Pointer transformFComp =
            CompositeTransformType::New();

//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName(filename1);
//  reader->Update();
//  ImagePointerType fixedImage = reader->GetOutput();
//
//  ReaderType::Pointer reader2 = ReaderType::New();
//  reader2->SetFileName(filename2);
//  reader2->Update();
//  ImagePointerType movingImage = reader2->GetOutput();

    const unsigned int imageSize = 6;
//  const unsigned int imageDimensionality = 2;

    ImageType::SizeType size;
    size.Fill(imageSize);
    ImageType::IndexType index;
    index.Fill(0);
    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index);
    ImageType::SpacingType spacing;
    spacing.Fill(1.0);
    ImageType::PointType origin;
    origin.Fill(0);
    ImageType::DirectionType direction;
    direction.SetIdentity();

    /* Create simple test images. */
    ImageType::Pointer fixedImage = ImageType::New();
    fixedImage->SetRegions(region);
    fixedImage->SetSpacing(spacing);
    fixedImage->SetOrigin(origin);
    fixedImage->SetDirection(direction);
    fixedImage->Allocate();

    ImageType::Pointer movingImage = ImageType::New();
    movingImage->SetRegions(region);
    movingImage->SetSpacing(spacing);
    movingImage->SetOrigin(origin);
    movingImage->SetDirection(direction);
    movingImage->Allocate();

    /* Fill images */
    ImageRegionIterator<ImageType> itFixed(fixedImage, region);
    itFixed.GoToBegin();
    unsigned int count = 1;
    while (!itFixed.IsAtEnd()) {
        itFixed.Set(count * count);
        count++;
        ++itFixed;
    }
    ImageRegionIteratorWithIndex<ImageType> itMoving(movingImage, region);
    itMoving.GoToBegin();
    count = 1;
    while (!itMoving.IsAtEnd()) {
        itMoving.Set(count * count);
        count++;
        ++itMoving;
    }

    VectorType zero;
    float def_value = 2.5;
    def_value = -0.5;
    zero.Fill(def_value);
    FieldType::Pointer field = FieldType::New();
    field->SetRegions(fixedImage->GetLargestPossibleRegion());
    field->SetSpacing(fixedImage->GetSpacing());
    field->SetOrigin(fixedImage->GetOrigin());
    field->SetDirection(fixedImage->GetDirection());
    field->Allocate();
    field->FillBuffer(zero);

    FieldType::Pointer fieldInv = FieldType::New();
//  def_value = 0;
    zero.Fill(def_value * (-1.0));
    fieldInv->SetRegions(fixedImage->GetLargestPossibleRegion());
    fieldInv->SetSpacing(fixedImage->GetSpacing());
    fieldInv->SetOrigin(fixedImage->GetOrigin());
    fieldInv->SetDirection(fixedImage->GetDirection());
    fieldInv->Allocate();
    fieldInv->FillBuffer(zero);

    zero.Fill(def_value * (1.0));
    transformMtranslation->Translate(zero);
    zero.Fill(def_value * (1.0));
    transformMtranslation2->Translate(zero);

    transformMdeformation->SetDeformationField(field);
    transformMdeformation->SetInverseDeformationField(fieldInv);

    transformMComp->AddTransform(transformMtranslation);
//  transformMComp->AddTransform(transformMdeformation);
    transformFComp->AddTransform(transformFId);

    typedef itk::ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<ImageType, ImageType> MetricType;
//   typedef DemonsImageToImageObjectMetric< ImageType, ImageType, ImageType > MetricType;

    typedef MetricType::Pointer MetricTypePointer;
    MetricTypePointer metric = MetricType::New();

    itk::Size<ImageDimension> neighborhood_radius;
    neighborhood_radius.Fill(1);

    metric->SetRadius(neighborhood_radius);

    metric->SetNumberOfSampling(1000);


    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(movingImage);
//  metric->SetFixedTransform(transformFComp);
//  metric->SetMovingTransform(transformMComp);

//  transformMtranslation

    metric->SetFixedTransform(transformFId);
//    metric->SetMovingTransform(transformMdeformation);
  metric->SetMovingTransform(transformMtranslation2);

    std::cout << "fixedImage:" << std::endl;
    PrintImage(fixedImage);

    std::cout << "movingImage:" << std::endl;
    PrintImage(movingImage);

//  std::cout << "fixed transform field:" << std::endl;
//  PrintImage(field);

//    std::cout << "moving transform field:" << std::endl;
//    PrintImage(field);

    // metric->SetMovingTransform( movingTransform );

    /* Initialize. */
    try {
        std::cout << "Calling Initialize..." << std::endl;
        metric->Initialize();
    } catch (ExceptionObject & exc) {
        std::cout << "Caught unexpected exception during Initialize: " << exc;
        return EXIT_FAILURE;
    }

    // Evaluate
    MetricType::MeasureType valueReturn;
    MetricType::DerivativeType derivativeReturn;
    try {
        std::cout << "Calling GetValueAndDerivative..." << std::endl;
        metric->GetValueAndDerivative(valueReturn, derivativeReturn);
    } catch (ExceptionObject & exc) {
        std::cout
                << "Caught unexpected exception during GetValueAndDerivative: "
                << exc;
        return EXIT_FAILURE;
    }

    std::cout << "Test passed." << std::endl;

    std::cout << "transformMdeformation parameters" << std::endl;

    std::cout << transformMdeformation->GetParameters() << std::endl;
    PrintImage(transformMdeformation->GetDeformationField());

//    std::cout << "transformMtranslation2 parameters" << std::endl;
//    std::cout << transformMtranslation2->GetParameters() << std::endl;

    std::cout << "derivative of moving transform:" << std::endl;
    std::cout << derivativeReturn << std::endl;

    std::cout << std::endl << "derivative of moving transform as a field:" << std::endl;
    PrintDerivativeAsVectorImage(fixedImage, derivativeReturn, ImageDimension);

//  PrintImage(fieldInv);

//  objectMetric->SetVirtualDomainSize(fixed_image->GetRequestedRegion().GetSize());
//  objectMetric->SetVirtualDomainIndex(fixed_image->GetRequestedRegion().GetIndex());
//  objectMetric->SetVirtualDomainSpacing(fixed_image->GetSpacing());
//  objectMetric->SetVirtualDomainOrigin(fixed_image->GetOrigin());
//  objectMetric->SetVirtualDomainDirection(fixed_image->GetDirection());
//
//


    return EXIT_SUCCESS;

}

//
//
//int main(int argc, char * argv[])
//{
//
//  if (argc <= 3) {
//    std::cout << "Args: 2D_fixed_image 2D_moving_image radius_in_one_dim number_of_threads" << std::endl;
//    exit(-1);
//  }
//
//  const char * filename1 = argv[1];
//  const char * filename2 = argv[2];
//  int radius_in_one_dim = atoi(argv[3]);
//  int number_of_threads = atoi(argv[4]);
//
//
//
//
//
//  const int ImageDimension = 2;
//  typedef itk::Image< unsigned char, ImageDimension > ImageType;
//  typedef ImageType::Pointer ImagePointerType;
//  typedef ImageType::RegionType RegionType;
//
//  typedef itk::Vector<float, ImageDimension> VectorType;
//  typedef itk::Image<VectorType, ImageDimension> VectorImageType;
//
//
//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName(filename1);
//  reader->Update();
//  ImagePointerType fixed_image = reader->GetOutput();
//
//  ReaderType::Pointer reader2 = ReaderType::New();
//  reader2->SetFileName(filename2);
//  reader2->Update();
//  ImagePointerType moving_image = reader2->GetOutput();
//
//  VectorImageType::Pointer field = VectorImageType::New();
//  field->SetRegions(fixed_image->GetLargestPossibleRegion());
//  field->Allocate();
//
//  VectorImageType::Pointer fieldInv = VectorImageType::New();
//  fieldInv->SetRegions(fixed_image->GetLargestPossibleRegion());
//  fieldInv->Allocate();
//
//  typedef itk::NeighborhoodNormalizedCrossCorrelationImageToImageMetric<ImageType, ImageType> ObjectMetricType;
//  typedef itk::MetricThreadedHolder<ObjectMetricType, VectorImageType> MetricThreadedHolderType;
//  typedef itk::ImageToData<ImageDimension, MetricThreadedHolderType> MetricThreaderType;
//
//
//  itk::Size<ImageDimension> neighborhood_radius;
//  neighborhood_radius.Fill(radius_in_one_dim);
//
//  // pseudo code
//  ObjectMetricType::Pointer objectMetric = ObjectMetricType::New();
//  MetricThreadedHolderType metricHolder;
//  MetricThreaderType::Pointer metricThreader = MetricThreaderType::New();
//
//
//  objectMetric->SetFixedImage(fixed_image);
//  objectMetric->SetMovingImage(moving_image);
//  objectMetric->SetRadius(neighborhood_radius);
//  objectMetric->InitializeGradientCalculator();
//
////  metricHolder.SetMetricFunction(objectMetric);
//  metricHolder.metric = objectMetric;
//  metricHolder.updateField = field;
//  metricHolder.updateFieldInv = fieldInv;
//  metricHolder.measure_per_thread.resize(number_of_threads);
//
//  ImageType::RegionType inboundary_region = fixed_image->GetLargestPossibleRegion();
//
//  metricThreader->SetNumberOfThreads(number_of_threads);
//  metricThreader->m_OverallRegion = inboundary_region ;
//  metricThreader->m_Holder = &metricHolder;
//  metricThreader->ThreadedGenerateData = MetricThreadedHolderType::ComputeMetricValueInRegionOnTheFlyThreaded;
//
//  metricThreader->GenerateData();
//
//  float energy = static_cast<float> (metricHolder.AccumulateMeasuresFromAllThreads());
//
//  std::cout << "cross correlation = " << energy << std::endl;
//
//  return 1;
//
//}
