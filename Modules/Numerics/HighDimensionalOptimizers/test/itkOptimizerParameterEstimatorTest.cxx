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
#include "itkTranslationTransform.h"

using namespace itk;

/* Create a simple metric to use for testing here. */
template< class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage >
class ITK_EXPORT OptimizerParameterEstimatorTestMetric:
  public itk::ObjectToObjectMetric
{
public:
  /** Standard class typedefs. */
  typedef OptimizerParameterEstimatorTestMetric                   Self;
  typedef itk::ObjectToObjectMetric                               Superclass;
  typedef itk::SmartPointer< Self >                               Pointer;
  typedef itk::SmartPointer< const Self >                         ConstPointer;

  typedef typename Superclass::MeasureType          MeasureType;
  typedef typename Superclass::DerivativeType       DerivativeType;
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  itkTypeMacro(OptimizerParameterEstimatorTestMetric, ObjectToObjectMetric);

  itkNewMacro(Self);

  // Pure virtual functions that all Metrics must provide
  unsigned int GetNumberOfParameters() const { return 5; }

  //using Superclass::GetValue;
  MeasureType GetValue()
    {
    return 1.0;
    }

  //using Superclass::GetValueAndDerivative;
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
  itkStaticConstMacro(FixedImageDimension, IndexValueType,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, IndexValueType,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, IndexValueType,
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

  OptimizerParameterEstimatorTestMetric() {}
  ~OptimizerParameterEstimatorTestMetric() {}

};

/**
 */
int itkOptimizerParameterEstimatorTest(int , char* [])
{

  // Image begins
  const IndexValueType  ImageDimension = 2;
  typedef double        PixelType;

  // Image Types
  typedef itk::Image<PixelType,ImageDimension>           FixedImageType;
  typedef itk::Image<PixelType,ImageDimension>           MovingImageType;
  typedef itk::Image<PixelType,ImageDimension>           VirtualImageType;

  FixedImageType::Pointer  fixedImage  = FixedImageType::New();
  MovingImageType::Pointer movingImage = MovingImageType::New();
  VirtualImageType::Pointer virtualImage = fixedImage;

  MovingImageType::SizeType    size;
  size.Fill(100);

  movingImage->SetRegions( size );
  fixedImage->SetRegions( size );
  // Image done

  // Transform begins
  typedef AffineTransform<double, ImageDimension>      MovingTransformType;
  MovingTransformType::Pointer movingTransform =  MovingTransformType::New();
  movingTransform->SetIdentity();

  typedef TranslationTransform<double, ImageDimension> FixedTransformType;
  FixedTransformType::Pointer fixedTransform =    FixedTransformType::New();
  fixedTransform->SetIdentity();
  // Transform done

  // Metric begins
  typedef OptimizerParameterEstimatorTestMetric<FixedImageType, MovingImageType>   MetricType;
  MetricType::Pointer metric = MetricType::New();

  metric->SetVirtualDomainImage( virtualImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );

  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );
  // Metric done

  // Testing OptimizerParameterEstimator

  typedef itk::OptimizerParameterEstimator< MetricType > OptimizerParameterEstimatorType;
  OptimizerParameterEstimatorType::Pointer parameterEstimator = OptimizerParameterEstimatorType::New();

  parameterEstimator->SetMetric(metric);
  parameterEstimator->Print( std::cout );

  // Scales for the moving transform
  parameterEstimator->SetTransformForward(true); //by default
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromShift); //by default
  OptimizerParameterEstimatorType::ScalesType movingScales(movingTransform->GetNumberOfParameters());

  parameterEstimator->EstimateScales(movingScales);
  std::cout << "Shift scales for the affine transform = " << movingScales << std::endl;

  // Check the correctness
  OptimizerParameterEstimatorType::ScalesType theoreticalMovingScales(movingTransform->GetNumberOfParameters());
  VirtualImageType::PointType upperPoint;
  virtualImage->TransformIndexToPhysicalPoint(virtualImage->GetLargestPossibleRegion().GetUpperIndex(), upperPoint);

  IndexValueType param = 0;
  for (IndexValueType row = 0; row < ImageDimension; row++)
    {
    for (IndexValueType col = 0; col < ImageDimension; col++)
      {
      theoreticalMovingScales[param++] = upperPoint[col] * upperPoint[col];
      }
    }
  for (IndexValueType row = 0; row < ImageDimension; row++)
    {
    theoreticalMovingScales[param++] = 1;
    }

  bool affinePass = true;
  for (IndexValueType p = 0; p < theoreticalMovingScales.GetSize(); p++)
    {
    if (movingScales[p] != theoreticalMovingScales[p])
      {
      affinePass = false;
      break;
      }
    }
  bool nonUniformForAffine = false;
  for (IndexValueType p = 1; p < movingScales.GetSize(); p++)
    {
    if (movingScales[p] != movingScales[0])
      {
      nonUniformForAffine = true;
      break;
      }
    }

  // Check done

  // Scales for the fixed transform
  parameterEstimator->SetTransformForward(false);
  OptimizerParameterEstimatorType::ScalesType fixedScales(fixedTransform->GetNumberOfParameters());

  parameterEstimator->EstimateScales(fixedScales);
  std::cout << "Shift scales for the translation transform = " << fixedScales << std::endl;

  // Check the correctness
  OptimizerParameterEstimatorType::ScalesType theoreticalFixedScales(fixedTransform->GetNumberOfParameters());
  theoreticalFixedScales.Fill(1.0);

  bool translationPass = true;
  for (IndexValueType p = 0; p < theoreticalFixedScales.GetSize(); p++)
    {
    if (fixedScales[p] != theoreticalFixedScales[p])
      {
      translationPass = false;
      break;
      }
    }
  bool uniformForTranslation = true;
  for (IndexValueType p = 1; p < fixedScales.GetSize(); p++)
    {
    if (fixedScales[p] != fixedScales[0])
      {
      uniformForTranslation = false;
      break;
      }
    }
  // Check done

  // Scales for the affine transform from transform jacobians
  parameterEstimator->SetTransformForward(true); //by default
  parameterEstimator->SetScaleStrategy(OptimizerParameterEstimatorType::ScalesFromJacobian);
  OptimizerParameterEstimatorType::ScalesType jacobianScales(movingTransform->GetNumberOfParameters());

  parameterEstimator->EstimateScales(jacobianScales);
  std::cout << "Jacobian scales for the affine transform = " << jacobianScales << std::endl;

  // Check the correctness
  OptimizerParameterEstimatorType::ScalesType theoreticalJacobianScales(movingTransform->GetNumberOfParameters());

  param = 0;
  for (IndexValueType row = 0; row < ImageDimension; row++)
    {
    for (IndexValueType col = 0; col < ImageDimension; col++)
      {
      //average of squares of consecutive integers [0,1,...,n]
      // = n*(n+1)*(2n+1)/6 / (n+1) = n*(2n+1)/6
      theoreticalJacobianScales[param++] = (upperPoint[col] * (2*upperPoint[col]+1)) / 6.0;
      }
    }
  for (IndexValueType row = 0; row < ImageDimension; row++)
    {
    theoreticalJacobianScales[param++] = 1;
    }

  bool jacobianPass = true;
  for (IndexValueType p = 0; p < jacobianScales.GetSize(); p++)
    {
    //due to random sampling, it is not exactly equal
    if (vcl_abs((jacobianScales[p] - theoreticalJacobianScales[p])/theoreticalJacobianScales[p]) > 0.2 )
      {
      jacobianPass = false;
      break;
      }
    }
  bool nonUniformForJacobian = false;
  for (IndexValueType p = 1; p < jacobianScales.GetSize(); p++)
    {
    if (jacobianScales[p] != jacobianScales[0])
      {
      nonUniformForJacobian = true;
      break;
      }
    }
  // Check done

  // Testing OptimizerParameterEstimator done
  std::cout << std::endl;

  if (!affinePass)
    {
    std::cout << "Failed: the shift scales for the affine transform are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the shift scales for the affine transform are correct." << std::endl;
    }

  if (!translationPass)
    {
    std::cout << "Failed: the shift scales for the translation transform are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the shift scales for the translation transform are correct." << std::endl;
    }

  if (!jacobianPass)
    {
    std::cout << "Failed: the jacobian scales for the affine transform are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the jacobian scales for the affine transform are correct." << std::endl;
    }

  if (!uniformForTranslation)
    {
    std::cout << "Error: the shift scales for a translation transform are not equal for all parameters." << std::endl;
    }
  if (!nonUniformForAffine)
    {
    std::cout << "Error: the shift scales for an affine transform are equal for all parameters." << std::endl;
    }
  if (!nonUniformForJacobian)
    {
    std::cout << "Error: the jacobian scales for an affine transform are equal for all parameters." << std::endl;
    }

  if (affinePass && translationPass && jacobianPass
    && nonUniformForAffine && uniformForTranslation && nonUniformForJacobian)
    {
    std::cout << "Test passed" << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "Test failed" << std::endl;
    return EXIT_FAILURE;
    }
}
