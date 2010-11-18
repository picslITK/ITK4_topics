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

#include "itkObjectToObjectMetric.h"

int itkObjectToObjectMetricTest(int argc,char* argv[])
{
  typedef itk::Image< unsigned char, 3 > ImageType;
  typedef itk::ObjectToObjectMetric< ImageType, ImageType> objectMetricType;
  //objectMetricType::Pointer objectMetric = objectMetricType::New();
  //objectMetric->SetNumberOfThreads(1);
  return EXIT_SUCCESS;
}
