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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkOptimizerParameterEstimator.h"
#include "itkObjectToObjectMetric.h"

#include "itkImage.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"

using namespace itk;

namespace{

/* Create a simple metric to use for testing here. */
template< class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage >
class ITK_EXPORT ObjectToObjectMetricProxy:
  public itk::ObjectToObjectMetric
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectMetricProxy                               Self;
  typedef itk::ObjectToObjectMetric                               Superclass;
  typedef itk::SmartPointer< Self >                               Pointer;
  typedef itk::SmartPointer< const Self >                         ConstPointer;

  typedef typename Superclass::MeasureType          MeasureType;
  typedef typename Superclass::DerivativeType       DerivativeType;
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  itkTypeMacro(ObjectToObjectMetricProxy, ObjectToObjectMetric);

  itkNewMacro(Self);

  // Pure virtual functions that all Metrics must provide
  unsigned int GetNumberOfParameters() const { return 5; }
  MeasureType GetValue()
    {
    return 1.0;
    }
  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative )
    {
    value = 1.0;
    derivative.Fill(0.0);
    }

  unsigned int GetNumberOfLocalParameters() const
  { return 0; }

  bool HasLocalSupport() const
  { return false; }

  void UpdateTransformParameters( DerivativeType &, ParametersValueType ) {}

  const ParametersType & GetParameters() const
  { return m_Parameters; }

  void Initialize(void) throw ( itk::ExceptionObject ) {}

  void PrintSelf(std::ostream& os, itk::Indent indent) const
  { Superclass::PrintSelf( os, indent ); }

  ParametersType  m_Parameters;

  // Image related types
  typedef TFixedImage                             FixedImageType;
  typedef TMovingImage                            MovingImageType;
  typedef TVirtualImage                           VirtualImageType;

  typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
  typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;
  typedef typename VirtualImageType::Pointer      VirtualImagePointer;

  /* Set/get images */
  /** Connect the Fixed Image.  */
  itkSetConstObjectMacro(FixedImage, FixedImageType);
  /** Get the Fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);
  /** Connect the Moving Image.  */
  itkSetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the Moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);
  /** Set all virtual domain image */
  itkSetObjectMacro(VirtualDomainImage, VirtualImageType);
  /** Get the virtual domain image */
  itkGetObjectMacro(VirtualDomainImage, VirtualImageType);

  /* Image dimension accessors */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, unsigned int,
      ::itk::GetImageDimension<VirtualImageType>::ImageDimension);

  /**  Type of the Transform Base classes */
  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )>  MovingTransformType;

  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )>  FixedTransformType;

  typedef typename FixedTransformType::Pointer        FixedTransformPointer;
  typedef typename MovingTransformType::Pointer       MovingTransformPointer;

  /** Connect the fixed transform. */
  itkSetObjectMacro(FixedTransform, FixedTransformType);
  /** Get a pointer to the fixed transform.  */
  itkGetConstObjectMacro(FixedTransform, FixedTransformType);
  /** Connect the moving transform. */
  itkSetObjectMacro(MovingTransform, MovingTransformType);
  /** Get a pointer to the moving transform.  */
  itkGetConstObjectMacro(MovingTransform, MovingTransformType);

private:

  FixedImageConstPointer  m_FixedImage;
  MovingImageConstPointer m_MovingImage;
  VirtualImagePointer     m_VirtualDomainImage;

  FixedTransformPointer   m_FixedTransform;
  MovingTransformPointer  m_MovingTransform;

  ObjectToObjectMetricProxy() {}
  ~ObjectToObjectMetricProxy() {}
};

/* global defines */
const int ImageDimension = 2;
typedef Image<double, ImageDimension>                    ImageType;
typedef ImageType::Pointer                               ImagePointerType;

}//namespace

/**
 */
int itkOptimizerParameterEstimatorTest(int , char* [])
{

  // Image begins
  const unsigned int  Dimension = 2;
  typedef double      PixelType;

  // Fixed Image Type
  typedef itk::Image<PixelType,Dimension>           FixedImageType;
  // Moving Image Type
  typedef itk::Image<PixelType,Dimension>           MovingImageType;

  MovingImageType::Pointer movingImage = MovingImageType::New();
  FixedImageType::Pointer  fixedImage  = FixedImageType::New();

  MovingImageType::SizeType    size;
  size.Fill(100);

  movingImage->SetRegions( size );
  fixedImage->SetRegions( size );
  // Image done

  // Transform begins
  typedef AffineTransform<double, Dimension>      MovingTransformType;
  typedef MovingTransformType::ParametersType     MovingParametersType;
  MovingTransformType::Pointer movingTransform =  MovingTransformType::New();
  movingTransform->SetIdentity();

  typedef IdentityTransform<double, Dimension>    FixedTransformType;
  FixedTransformType::Pointer fixedTransform =    FixedTransformType::New();
  fixedTransform->SetIdentity();
  // Transform done

  // Metric begins
  typedef ObjectToObjectMetricProxy<ImageType,ImageType>   MetricType;
  MetricType::Pointer metric = MetricType::New();

  metric->SetVirtualDomainImage( fixedImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );
  // Metric done

  // Testing OptimizerParameterEstimator

  typedef itk::OptimizerParameterEstimator< MetricType > OptimizerParameterEstimatorType;
  OptimizerParameterEstimatorType::Pointer parameterEstimator = OptimizerParameterEstimatorType::New();

  parameterEstimator->SetMetric(metric);
  //parameterEstimator->Print( std::cout );

  parameterEstimator->SetTransformForward(true); //by default
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromShift); //by default
  //parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromJacobian);
  OptimizerParameterEstimatorType::ScalesType movingScales(metric->GetMovingTransform()->GetNumberOfParameters());

  parameterEstimator->EstimateScales(movingScales);
  std::cout << "Scales for moving transform parameters = " << movingScales << std::endl;

  parameterEstimator->SetTransformForward(false);
  OptimizerParameterEstimatorType::ScalesType fixedScales(metric->GetFixedTransform()->GetNumberOfParameters());

  parameterEstimator->EstimateScales(fixedScales);
  std::cout << "Scales for fixed transform parameters = " << fixedScales << std::endl;

  // Testing OptimizerParameterEstimator done

  std::cout << std::endl << "Test passed" << std::endl;

  return EXIT_SUCCESS;
}
