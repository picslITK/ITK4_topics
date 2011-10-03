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
#ifndef __itkRegistrationParameterScalesEstimator_h
#define __itkRegistrationParameterScalesEstimator_h

#include "itkTransform.h"
#include "itkOptimizerParameterScalesEstimator.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/** \class RegistrationParameterScalesEstimator
 *  \brief Implements a registration helper class for estimating scales of
 * transform parameters.
 *
 * Its input includes the fixed/moving images and transform objects,
 * which are obtained from the metric object.
 *
 * This class implements some common methods as building blocks called by
 * subclasses with various estimation strategies. One of these methods is
 * SampleImageDomain, which provides various choices of sampling the image
 * domain.
 *
 * \ingroup ITKHighDimensionalOptimizers
 */
template < class TMetric >
class ITK_EXPORT RegistrationParameterScalesEstimator : public OptimizerParameterScalesEstimator
{
public:
  /** Standard class typedefs. */
  typedef RegistrationParameterScalesEstimator  Self;
  typedef OptimizerParameterScalesEstimator     Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( RegistrationParameterScalesEstimator, OptimizerParameterScalesEstimator );

  /** Type of scales */
  typedef typename Superclass::ScalesType           ScalesType;
  /** Type of paramters of the optimizer */
  typedef typename Superclass::ParametersType       ParametersType;
  /** Type of float */
  typedef typename Superclass::FloatType            FloatType;

  typedef TMetric                                   MetricType;
  typedef typename MetricType::Pointer              MetricPointer;

  /** Type of the transform to initialize */
  typedef typename MetricType::FixedTransformType   FixedTransformType;
  typedef typename FixedTransformType::ConstPointer FixedTransformPointer;

  typedef typename MetricType::MovingTransformType  MovingTransformType;
  typedef typename MovingTransformType::ConstPointer
                                                    MovingTransformPointer;

  /** Image Types to use in the initialization of the transform */
  typedef typename TMetric::FixedImageType          FixedImageType;
  typedef typename TMetric::MovingImageType         MovingImageType;
  typedef typename TMetric::VirtualImageType        VirtualImageType;

  typedef typename FixedImageType::ConstPointer     FixedImagePointer;
  typedef typename MovingImageType::ConstPointer    MovingImagePointer;
  typedef typename VirtualImageType::ConstPointer   VirtualImagePointer;

  /* Image dimension accessors */
  itkStaticConstMacro(FixedImageDimension, SizeValueType,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, SizeValueType,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, SizeValueType,
      ::itk::GetImageDimension<VirtualImageType>::ImageDimension);

  typedef typename VirtualImageType::RegionType     VirtualRegionType;
  typedef typename VirtualImageType::SizeType       VirtualSizeType;
  typedef typename VirtualImageType::PointType      VirtualPointType;
  typedef typename VirtualImageType::IndexType      VirtualIndexType;

  typedef typename FixedImageType::PointType        FixedPointType;
  typedef typename FixedImageType::IndexType        FixedIndexType;
  typedef typename FixedImageType::PointValueType   FixedPointValueType;
  typedef typename itk::ContinuousIndex< FixedPointValueType,
          FixedImageType::ImageDimension >          FixedContinuousIndexType;

  typedef typename MovingImageType::PointType       MovingPointType;
  typedef typename MovingImageType::IndexType       MovingIndexType;
  typedef typename MovingImageType::PointValueType  MovingPointValueType;
  typedef typename itk::ContinuousIndex< MovingPointValueType,
          MovingImageType::ImageDimension >         MovingContinuousIndexType;

  /** The strategies to sample physical points in the virtual domain. */
  typedef enum { FullDomainSampling, CornerSampling, RandomSampling}
          SamplingStrategyType;

  typedef std::vector<VirtualPointType>             ImageSampleContainerType;

  /** Type of Jacobian of transform */
  typedef typename TMetric::FixedTransformJacobianType       FixedJacobianType;
  typedef typename TMetric::MovingTransformJacobianType      MovingJacobianType;

  /** Set the sampling strategy */
  itkSetMacro(SamplingStrategy, SamplingStrategyType);

  /** SetMetric sets the metric used in the estimation process.
   *  The images and transforms from the metric will be used for estimation.
   */
  itkSetObjectMacro(Metric, MetricType);

  /** m_TransformForward specifies which transform scales to be estimated.
   * m_TransformForward = true for the moving transform parameters.
   * m_TransformForward = false for the fixed transform parameters.
   */
  itkSetMacro(TransformForward, bool);
  itkGetConstMacro(TransformForward, bool);

  /* Set and get the number of image samples. */
  itkSetMacro(NumberOfRandomSamples, SizeValueType);

  /** Estimate parameter scales */
  virtual void EstimateScales(ScalesType &scales) = 0;

protected:
  RegistrationParameterScalesEstimator();
  ~RegistrationParameterScalesEstimator(){};

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  /** Check and set the images and transforms from the metric. */
  bool CheckAndSetInputs();

  /** Transform a point to its continous index. */
  template< class TContinuousIndexType > void TransformPointToContinuousIndex(
                              const VirtualPointType &point,
                              TContinuousIndexType &mappedIndex);

  /** Compute the transform Jacobian at a physical point. */
  template< class TJacobianType > void ComputeSquaredJacobianNorms(
                              const VirtualPointType  & p,
                              ParametersType & squareNorms);

  /** Sample the virtual image domain and store the physical points in m_ImageSamples. */
  virtual void SampleImageDomain();

  /** Sample the virutal domain with all pixels. */
  void SampleImageDomainFully();

  /** Sample the virutal domain with corners. */
  void SampleImageDomainWithCorners();

  /** Sample the virutal domain randomly in a uniform distribution. */
  void SampleImageDomainRandomly();

  /** Get the transform in use. */
  const TransformBase *GetTransform();

  /** Get the dimension of the target image transformed to. */
  SizeValueType GetImageDimension();

  /** Get the fixed Image. */
  itkGetConstObjectMacro(FixedImage,  FixedImageType);
  /** Get the moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the virtual Image. */
  itkGetConstObjectMacro(VirtualImage, VirtualImageType);

  /** Get the fixed transform. */
  itkGetConstObjectMacro(FixedTransform,  FixedTransformType);
  /** Get the moving transform. */
  itkGetConstObjectMacro(MovingTransform, MovingTransformType);

  // the image samples in the virtual image domain
  ImageSampleContainerType      m_ImageSamples;

  // the number of image samples in the virtual image domain
  SizeValueType                 m_NumberOfRandomSamples;

  static const SizeValueType    SizeOfSmallDomain = 1000;

private:
  RegistrationParameterScalesEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // the metric object
  MetricPointer                 m_Metric;

  // the transform objects
  FixedTransformPointer         m_FixedTransform;
  MovingTransformPointer        m_MovingTransform;

  // the fixed images
  FixedImagePointer             m_FixedImage;
  // the moving images
  MovingImagePointer            m_MovingImage;
  // the virtual image for symmetric registration
  VirtualImagePointer           m_VirtualImage;

  /** m_TransformForward specifies which transform scales to be estimated.
   * m_TransformForward = true for the moving transform parameters.
   * m_TransformForward = false for the fixed transform parameters.
   */
  bool m_TransformForward;

  // sampling stategy
  SamplingStrategyType          m_SamplingStrategy;

}; //class RegistrationParameterScalesEstimator


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRegistrationParameterScalesEstimator.hxx"
#endif

#endif /* __itkRegistrationParameterScalesEstimator_h */
