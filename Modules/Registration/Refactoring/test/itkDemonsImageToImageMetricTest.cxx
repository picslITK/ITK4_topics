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
#include "itkObjectToObjectMetric.h"
#include "itkObjectToObjectThreadedMetricOptimizer.h"
#include "itkDemonsImageToImageMetric.h"
#include "itkIdentityTransform.h"
#include "itkDeformationFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"
#include "itkImageToData.h"


int itkDemonsImageToImageMetricTest(int argc, char * argv[])
{

  if (argc <= 3) {
    std::cout << "Args: 2D_fixed_image 2D_moving_image number_of_threads" << std::endl;
    exit(-1);
  }

  const char * filename1 = argv[1];
  const char * filename2 = argv[2];
  int number_of_threads = atoi(argv[3]);

  const int ImageDimension = 2;
  //  typedef itk::VectorImage< float,ImageDimension > ImageType;
  typedef itk::Image<float, ImageDimension> ImageType;
  typedef ImageType::Pointer ImagePointerType;
  typedef ImageType::RegionType RegionType;

  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  typedef itk::Transform<double,ImageDimension>  TransformType;
  typedef itk::IdentityTransform<double,ImageDimension>  IdentityTransformType;
  typedef itk::CompositeTransform<double,ImageDimension>  CompositeTransformType;
  typedef itk::TranslationTransform<double,ImageDimension>  TranslationTransformType;
  typedef itk::DeformationFieldTransform<double,ImageDimension>  DeformationTransformType;
  typedef DeformationTransformType::DeformationFieldType FieldType;
  IdentityTransformType::Pointer transformFId = IdentityTransformType::New();
  IdentityTransformType::Pointer transformMId = IdentityTransformType::New();
  DeformationTransformType::Pointer transformMdeformation = DeformationTransformType::New();
  TranslationTransformType::Pointer transformMtranslation = TranslationTransformType::New();
  TranslationTransformType::Pointer transformMtranslation2 = TranslationTransformType::New();
  CompositeTransformType::Pointer transformMComp = CompositeTransformType::New();
  CompositeTransformType::Pointer transformFComp = CompositeTransformType::New();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename1);
  reader->Update();
  ImagePointerType fixed_image = reader->GetOutput();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName(filename2);
  reader2->Update();
  ImagePointerType moving_image = reader2->GetOutput();

  VectorType zero;
  float def_value=2.5;
  def_value=0.5;
  zero.Fill(def_value);
  FieldType::Pointer field = FieldType::New();
  field->SetRegions(fixed_image->GetLargestPossibleRegion());
  field->SetSpacing(fixed_image->GetSpacing());
  field->SetOrigin(fixed_image->GetOrigin());
  field->SetDirection(fixed_image->GetDirection());
  field->Allocate();
  field->FillBuffer(zero);

  FieldType::Pointer fieldInv = FieldType::New();
  zero.Fill(def_value*(-1.0));
  fieldInv->SetRegions(fixed_image->GetLargestPossibleRegion());
  fieldInv->SetSpacing(fixed_image->GetSpacing());
  fieldInv->SetOrigin(fixed_image->GetOrigin());
  fieldInv->SetDirection(fixed_image->GetDirection());
  fieldInv->Allocate();
  fieldInv->FillBuffer(zero);

  zero.Fill(def_value*(1.0));
  transformMtranslation->Translate(zero);
  zero.Fill(def_value*(1.0));
  transformMtranslation2->Translate(zero);

  transformMdeformation->SetDeformationField(field);
  transformMdeformation->SetInverseDeformationField(fieldInv);

  transformMComp->AddTransform(transformMtranslation);
//  transformMComp->AddTransform(transformMdeformation);
  transformFComp->AddTransform(transformFId);
  typedef itk::DemonsImageToImageMetric<ImageType,ImageType> ObjectMetricType;
  typedef ObjectMetricType::Pointer MetricTypePointer;
  MetricTypePointer objectMetric = ObjectMetricType::New();
  objectMetric->SetFixedImage(fixed_image);
  objectMetric->SetMovingImage(moving_image);
  objectMetric->SetFixedImageTransform(transformFComp);
  objectMetric->SetMovingImageTransform(transformMComp);
  objectMetric->SetVirtualDomainSize(fixed_image->GetRequestedRegion().GetSize());
  objectMetric->SetVirtualDomainIndex(fixed_image->GetRequestedRegion().GetIndex());
  objectMetric->SetVirtualDomainSpacing(fixed_image->GetSpacing());
  objectMetric->SetVirtualDomainOrigin(fixed_image->GetOrigin());
  objectMetric->SetVirtualDomainDirection(fixed_image->GetDirection());
  // Compute one iteration of the metric
  objectMetric->Initialize();

  /* The threader type might be definable w/in the optimizer ? */
  typedef itk::ObjectToObjectThreadedMetricOptimizer<ObjectMetricType>
                                                  MetricThreadedOptimizerType;
  itk::Size<ImageDimension> neighborhood_radius;
  neighborhood_radius.Fill(0);

  // pseudo code
  MetricThreadedOptimizerType::Pointer metricOptimizer =
    MetricThreadedOptimizerType::New();
  metricOptimizer->SetMetric( objectMetric );
  ImageType::RegionType inboundary_region = fixed_image->GetRequestedRegion();
  metricOptimizer->SetOverallRegion( inboundary_region );
  metricOptimizer->SetNumberOfThreads(number_of_threads);

  metricOptimizer->StartOptimization();

  float energy = static_cast<float> ( metricOptimizer->GetValue() );

  std::cout << "metric = " << energy << std::endl;
  return EXIT_FAILURE; //not yet a full test

}
