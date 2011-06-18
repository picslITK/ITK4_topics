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
#include "itkImage.h"
#include "itkVector.h"
#include "itkImageToImageObjectMetric.h"
#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"

using namespace itk;

namespace {

template<class TFixedImage,class TMovingImage,class TVirtualImage>
class TestDerivedMetric
  : public ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:
  /** Standard class typedefs. */
  typedef TestDerivedMetric                                   Self;
  typedef ImageToImageObjectMetric<TFixedImage, TMovingImage,
                                               TVirtualImage> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /** superclass types */
  typedef typename Superclass::MeasureType                    MeasureType;
  typedef typename Superclass::DerivativeType                 DerivativeType;

  /* Implement pure virtual methods */
  void Initialize() throw ( itk::ExceptionObject )
  {
    //Be sure to call base class initialize
    Superclass::Initialize();

    //Now do your own initialization here
  }

  MeasureType GetValue()
  {
    //TODO
    return 1.0;
  }

  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative)
  {
    //1) Do any pre-processing required for your metric

    //2) Call GetValueAndDerivativeMultiThreadedInitiate.
    //This will iterate of image region and call your
    // GetValueAndDerivativeProcessPoint method, see definition in
    // base.

    //3) Optionally call GetValueAndDerivativeMultiThreadedPostProcess for
    // default post-processing, which sums up results from each thread,
    // and optionally averages them.
    //Do your own post-processing as needed.

    //That's it. Easy as 1, 2, 3.
  }

protected:
  TestDerivedMetric(){};
  virtual ~TestDerivedMetric() {}
  void PrintSelf(std::ostream& os, Indent indent) const {}

private:
  //purposely not implemented
  TestDerivedMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);
};
}//namespace

////////////////////////////////////////////////////////////

int itkImageToImageObjectMetricTest(int argc, char * argv[])
{
  typedef Image< double, 2 > ImageType;
  typedef TestDerivedMetric<ImageType,ImageType,ImageType> TestMetricType;
  TestMetricType::Pointer metric = TestMetricType::New();

  return EXIT_SUCCESS;
}
