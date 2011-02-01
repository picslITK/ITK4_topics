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

// #include "itkObjectToObjectMetric.h"

#include "itkDemonsImageToImageMetric.h"
#include "itkIdentityTransform.h"
#include "itkDeformationFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"
#include "itkImageToData.h"


int main(int argc, char * argv[])
{

  if (argc <= 3) {
    std::cout << "Args: 2D_fixed_image 2D_moving_image number_of_threads" << std::endl;
    exit(-1);
  }

  const char * filename1 = argv[1];
  const char * filename2 = argv[2];
  //  int number_of_threads = atoi(argv[3]);

  const int ImageDimension = 2;
  typedef itk::Image< unsigned char, ImageDimension > ImageType;
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
  IdentityTransformType::Pointer transformF = IdentityTransformType::New();
  DeformationTransformType::Pointer transformM1 = DeformationTransformType::New();
  TranslationTransformType::Pointer transformM2 = TranslationTransformType::New();
  TranslationTransformType::Pointer transformM3 = TranslationTransformType::New();
  CompositeTransformType::Pointer transformM = CompositeTransformType::New();

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
  float def_value=1.0;
  zero.Fill(def_value);
  FieldType::Pointer field = FieldType::New();
  field->SetRegions(fixed_image->GetLargestPossibleRegion());
  field->Allocate();
  field->FillBuffer(zero);

  FieldType::Pointer fieldInv = FieldType::New();
  zero.Fill(def_value*(-1.0));
  fieldInv->SetRegions(fixed_image->GetLargestPossibleRegion());
  fieldInv->Allocate();
  fieldInv->FillBuffer(zero);

  zero.Fill(def_value*(1.0));
  transformM2->Translate(zero);
  zero.Fill(def_value*(1.0));
  transformM3->Translate(zero);

  transformM1->SetDeformationField(field);
  transformM1->SetInverseDeformationField(field);
  
  transformM->AddTransform(transformM2);
  transformM->AddTransform(transformM1);

  typedef itk::DemonsImageToImageMetric<ImageType, ImageType> ObjectMetricType;
  ObjectMetricType::Pointer objectMetric = ObjectMetricType::New();

  /** Set up the virtual reference space */
  objectMetric->SetVirtualDomainSpacing(fixed_image->GetSpacing());
  objectMetric->SetVirtualDomainSize(fixed_image->GetLargestPossibleRegion().GetSize());
  objectMetric->SetVirtualDomainOrigin(fixed_image->GetOrigin());
  objectMetric->SetVirtualDomainDirection(fixed_image->GetDirection());

  /** Define the input images and their transforms */
  objectMetric->SetFixedImage(fixed_image);
  objectMetric->SetMovingImage(moving_image);
  objectMetric->SetFixedImageTransform(transformF);
  objectMetric->SetMovingImageTransform(transformM);

  /** Compute one iteration of the metric */
  objectMetric->ComputeMetricAndDerivative();

  return 1;

}
