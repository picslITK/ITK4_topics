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
#include "itkRegistrationParameterScalesEstimator.h"
#include "itkObjectToObjectMetric.h"

#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkTranslationTransform.h"

/**
 *  \class RegistrationParameterScalesEstimatorTestMetric for test.
 *  Create a simple metric to use for testing here.
 */
template< class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage >
class ITK_EXPORT RegistrationParameterScalesEstimatorTestMetric:
  public itk::ObjectToObjectMetric
{
public:
  /** Standard class typedefs. */
  typedef RegistrationParameterScalesEstimatorTestMetric          Self;
  typedef itk::ObjectToObjectMetric                               Superclass;
  typedef itk::SmartPointer< Self >                               Pointer;
  typedef itk::SmartPointer< const Self >                         ConstPointer;

  typedef typename Superclass::MeasureType          MeasureType;
  typedef typename Superclass::DerivativeType       DerivativeType;
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  itkTypeMacro(RegistrationParameterScalesEstimatorTestMetric, ObjectToObjectMetric);

  itkNewMacro(Self);

  // Pure virtual functions that all Metrics must provide
  unsigned int GetNumberOfParameters() const { return 5; }

  MeasureType GetValue() const
    {
    return 1.0;
    }

  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
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
  itkStaticConstMacro(FixedImageDimension, itk::SizeValueType,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, itk::SizeValueType,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, itk::SizeValueType,
      ::itk::GetImageDimension<VirtualImageType>::ImageDimension);

  /**  Type of the Transform Base classes */
  typedef ::itk::Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )>  MovingTransformType;

  typedef ::itk::Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )>  FixedTransformType;

  typedef typename FixedTransformType::Pointer        FixedTransformPointer;
  typedef typename MovingTransformType::Pointer       MovingTransformPointer;

  typedef typename FixedTransformType::JacobianType   FixedTransformJacobianType;
  typedef typename MovingTransformType::JacobianType  MovingTransformJacobianType;

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

  RegistrationParameterScalesEstimatorTestMetric() {}
  ~RegistrationParameterScalesEstimatorTestMetric() {}

};

/**
 *  \class RegistrationParameterScalesEstimatorTest for test.
 *  Create a simple scales estimator class to use for testing here.
 */
template < class TMetric >
class ITK_EXPORT RegistrationParameterScalesEstimatorTest:
  public itk::RegistrationParameterScalesEstimator< TMetric >
{
public:
  /** Standard class typedefs. */
  typedef RegistrationParameterScalesEstimatorTest                    Self;
  typedef itk::RegistrationParameterScalesEstimator< TMetric >        Superclass;
  typedef itk::SmartPointer< Self >                                   Pointer;
  typedef itk::SmartPointer< const Self >                             ConstPointer;

  itkNewMacro(Self);

  itkTypeMacro(RegistrationParameterScalesEstimatorTest, RegistrationParameterScalesEstimator);

  /** Type of scales */
  typedef typename Superclass::ScalesType                ScalesType;
  /** Type of paramters of the optimizer */
  typedef typename Superclass::ParametersType            ParametersType;
  /** Type of float */
  typedef typename Superclass::FloatType                 FloatType;

  typedef typename Superclass::VirtualPointType          VirtualPointType;
  typedef typename Superclass::MovingTransformType       MovingTransformType;
  typedef typename Superclass::FixedTransformType        FixedTransformType;
  typedef typename Superclass::MovingJacobianType        MovingJacobianType;
  typedef typename Superclass::FixedJacobianType         FixedJacobianType;

  /** Estimate parameter scales with maximum squared norms of Jacobians. */
  virtual void EstimateScales(ScalesType &parameterScales)
    {
    this->CheckAndSetInputs();
    this->SampleImageDomain();

    itk::SizeValueType numPara = this->GetTransform()->GetNumberOfParameters();
    parameterScales.SetSize(numPara);

    ParametersType norms(numPara);

    itk::SizeValueType numSamples = this->m_ImageSamples.size();

    norms.Fill(0.0);
    parameterScales.Fill(1.0);

    // checking each sample point
    for (itk::SizeValueType c=0; c<numSamples; c++)
      {
      VirtualPointType point = this->m_ImageSamples[c];

      ParametersType squaredNorms(numPara);
      if (this->GetTransformForward())
        {
        this->template ComputeSquaredJacobianNorms<MovingJacobianType>( point, squaredNorms );
        }
      else
        {
        this->template ComputeSquaredJacobianNorms<FixedJacobianType>( point, squaredNorms );
        }
      for (itk::SizeValueType p=0; p<numPara; p++)
        {
        if (norms[p] < squaredNorms[p])
          {
          norms[p] = squaredNorms[p];
          }
        }
      } //for numSamples

    if (numSamples > 0)
      {
      for (itk::SizeValueType p=0; p<numPara; p++)
        {
        parameterScales[p] = norms[p];
        }
      }
    }

protected:
  RegistrationParameterScalesEstimatorTest(){};
  ~RegistrationParameterScalesEstimatorTest(){};

private:
  RegistrationParameterScalesEstimatorTest(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

/**
 */
int itkRegistrationParameterScalesEstimatorTest(int , char* [])
{

  // Image begins
  const itk::SizeValueType ImageDimension = 2;
  typedef double           PixelType;

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
  typedef itk::AffineTransform<double, ImageDimension>      MovingTransformType;
  MovingTransformType::Pointer movingTransform =  MovingTransformType::New();
  movingTransform->SetIdentity();

  typedef itk::TranslationTransform<double, ImageDimension> FixedTransformType;
  FixedTransformType::Pointer fixedTransform =    FixedTransformType::New();
  fixedTransform->SetIdentity();
  // Transform done

  // Metric begins
  typedef RegistrationParameterScalesEstimatorTestMetric
    <FixedImageType, MovingImageType> MetricType;
  MetricType::Pointer metric = MetricType::New();

  metric->SetVirtualDomainImage( virtualImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );

  metric->SetFixedTransform( fixedTransform );
  metric->SetMovingTransform( movingTransform );
  // Metric done

  // Scales for the affine transform from max squared norm of transform jacobians
  typedef RegistrationParameterScalesEstimatorTest< MetricType >
    RegistrationParameterScalesEstimatorTestType;
  RegistrationParameterScalesEstimatorTestType::Pointer jacobianScaleEstimator
    = RegistrationParameterScalesEstimatorTestType::New();

  jacobianScaleEstimator->SetMetric(metric);
  jacobianScaleEstimator->SetTransformForward(true);
  jacobianScaleEstimator->SetSamplingStrategy(
    RegistrationParameterScalesEstimatorTestType::CornerSampling);
  jacobianScaleEstimator->Print( std::cout );

  RegistrationParameterScalesEstimatorTestType::ScalesType jacobianScales(
    movingTransform->GetNumberOfParameters());
  jacobianScaleEstimator->EstimateScales(jacobianScales);
  std::cout << "Scales from max squared Jacobian norm for the affine transform = "
    << jacobianScales << std::endl;

  // Check the correctness
  RegistrationParameterScalesEstimatorTestType::ScalesType theoreticalJacobianScales(
    movingTransform->GetNumberOfParameters());
  VirtualImageType::PointType upperPoint;
  virtualImage->TransformIndexToPhysicalPoint(virtualImage->
    GetLargestPossibleRegion().GetUpperIndex(), upperPoint);

  itk::SizeValueType param = 0;
  for (itk::SizeValueType row = 0; row < ImageDimension; row++)
    {
    for (itk::SizeValueType col = 0; col < ImageDimension; col++)
      {
      // max squared jacobian norms
      theoreticalJacobianScales[param++] = upperPoint[col] * upperPoint[col];
      }
    }
  for (itk::SizeValueType row = 0; row < ImageDimension; row++)
    {
    theoreticalJacobianScales[param++] = 1;
    }

  bool jacobianPass = true;
  for (itk::SizeValueType p = 0; p < jacobianScales.GetSize(); p++)
    {
    if (jacobianScales[p] != theoreticalJacobianScales[p])
      {
      jacobianPass = false;
      break;
      }
    }
  bool nonUniformForJacobian = false;
  for (itk::SizeValueType p = 1; p < jacobianScales.GetSize(); p++)
    {
    if (jacobianScales[p] != jacobianScales[0])
      {
      nonUniformForJacobian = true;
      break;
      }
    }
  // Check done

  // Testing different sampling strategies
  jacobianScaleEstimator->SetSamplingStrategy(
    RegistrationParameterScalesEstimatorTestType::RandomSampling);
  jacobianScaleEstimator->SetNumberOfRandomSamples( 1000 );
  jacobianScaleEstimator->EstimateScales(jacobianScales);
  bool randomPass = true;
  for (itk::SizeValueType p = 0; p < jacobianScales.GetSize(); p++)
    {
    if (vcl_abs( (jacobianScales[p] - theoreticalJacobianScales[p])
      / theoreticalJacobianScales[p] ) > 0.3 )
      {
      randomPass = false;
      break;
      }
    }
  jacobianScaleEstimator->SetSamplingStrategy(
    RegistrationParameterScalesEstimatorTestType::FullDomainSampling);
  jacobianScaleEstimator->EstimateScales(jacobianScales);
  bool fullDomainPass = true;
  for (itk::SizeValueType p = 0; p < jacobianScales.GetSize(); p++)
    {
    if (jacobianScales[p] != theoreticalJacobianScales[p])
      {
      fullDomainPass = false;
      break;
      }
    }

  // Testing RegistrationParameterScalesEstimatorTest done
  std::cout << std::endl;

  if (!jacobianPass)
    {
    std::cout << "Failed: the jacobian scales for the affine transform are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the jacobian scales for the affine transform are correct." << std::endl;
    }

  if (!randomPass)
    {
    std::cout << "Failed: the jacobian scales with random sampling are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the jacobian scales with random sampling are correct." << std::endl;
    }

  if (!fullDomainPass)
    {
    std::cout << "Failed: the jacobian scales from checking the full domain are not correct." << std::endl;
    }
  else
    {
    std::cout << "Passed: the jacobian scales from checking the full domain are correct." << std::endl;
    }

  if (!nonUniformForJacobian)
    {
    std::cout << "Error: the jacobian scales for an affine transform are equal for all parameters." << std::endl;
    }

  if (jacobianPass && nonUniformForJacobian && randomPass && fullDomainPass)
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
