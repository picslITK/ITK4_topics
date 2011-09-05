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
#ifndef __itkOptimizerParameterEstimator_h
#define __itkOptimizerParameterEstimator_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkTransformBase.h"
#include "itkOptimizerParameterEstimatorBase.h"
#include "itkImageRandomConstIteratorWithIndex.h"

#include <iostream>

namespace itk
{

/** \class OptimizerParameterEstimator
 *  \brief Implements an optimizer helper class for estimating scales of
 * transform parameters and computing the maximum voxel shift from a
 * specific transform. The maximum voxel shift may be used in an optimizer
 * to decide the step size.
 *
 * Its input includes the fixed/moving images and transform objects,
 * which can be got from the metric object.
 *
 * We've implemented two strategies to estimate scales. To choose one,
 * please call SetScaleStrategy(option) where the argument is either
 * ScalesFromShift or ScalesFromJacobian.
 *
 * With ScalesFromShift, the scale of a parameter is estimated from the
 * maximum voxel shift produced from a unit variation of this parameter. The
 * maximization is done by checking the corners of the image domain.
 * Checking the corners is enough for affine transforms, but may be changed
 * for other types of transforms.
 *
 * With the second strategy ScalesFromJacobian, the scale of a parameter is
 * estimated from the averaged squared norm of the Jacobian w.r.t the
 * parameter. The averaging is done over a uniformly random sampling of the
 * image domain.
 *
 * \ingroup ITKHighDimensionalOptimizers
 */
template < class TMetric >
class ITK_EXPORT OptimizerParameterEstimator : public OptimizerParameterEstimatorBase
{
public:
  /** Standard class typedefs. */
  typedef OptimizerParameterEstimator           Self;
  typedef OptimizerParameterEstimatorBase       Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimizerParameterEstimator, OptimizerParameterEstimatorBase );

  typedef TMetric                                   MetricType;
  typedef typename MetricType::Pointer              MetricPointer;

  /** Type of the transform to initialize */
  typedef typename MetricType::FixedTransformType   FixedTransformType;
  typedef typename FixedTransformType::Pointer      FixedTransformPointer;

  typedef typename MetricType::MovingTransformType  MovingTransformType;
  typedef typename MovingTransformType::Pointer     MovingTransformPointer;

  /** Image Types to use in the initialization of the transform */
  typedef typename TMetric::FixedImageType          FixedImageType;
  typedef typename TMetric::MovingImageType         MovingImageType;
  typedef typename TMetric::VirtualImageType        VirtualImageType;

  typedef typename FixedImageType::ConstPointer     FixedImagePointer;
  typedef typename MovingImageType::ConstPointer    MovingImagePointer;
  typedef typename VirtualImageType::ConstPointer   VirtualImagePointer;

  /* Image dimension accessors */
  itkStaticConstMacro(FixedImageDimension, IndexValueType,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, IndexValueType,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, IndexValueType,
      ::itk::GetImageDimension<VirtualImageType>::ImageDimension);

  /** The stratigies to decide scales */
  typedef enum { ScalesFromShift, ScalesFromJacobian } ScaleStrategyType;
  /** Set the learning rate strategy */
  itkSetMacro(ScaleStrategy, ScaleStrategyType);

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

  /** SetMetric sets the metric used in the estimation process. SetMetric
   *  gets the images and transforms from the metric. Please make sure the metric
   *  has these members set when SetMetric(metric) is called.
   */
  virtual void SetMetric(MetricType *metric);

  /** Check if the metric, images and transforms are properly set */
  bool CheckInputs() const;

  /** Get the metric. */
  itkGetConstObjectMacro(Metric, MetricType);

  /** Set the fixed image used in the registration process */
  itkSetConstObjectMacro(FixedImage,  FixedImageType);
  /** Get the fixed Image. */
  itkGetConstObjectMacro(FixedImage,  FixedImageType);

  /** Set the moving image used in the registration process */
  itkSetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Set the virtual image used in the registration process */
  itkSetConstObjectMacro(VirtualImage, VirtualImageType);
  /** Get the virtual Image. */
  itkGetConstObjectMacro(VirtualImage, VirtualImageType);

  /** Set the fixed transform */
  itkSetObjectMacro(FixedTransform,  FixedTransformType);
  /** Get the fixed transform */
  itkGetConstObjectMacro(FixedTransform,  FixedTransformType);

  /** Set the moving transform */
  itkSetObjectMacro(MovingTransform, MovingTransformType);
  /** Get the moving transform */
  itkGetConstObjectMacro(MovingTransform, MovingTransformType);

  /** m_TransformForward specifies which transform scales to be estimated.
   * m_TransformForward = true for the moving transform parameters.
   * m_TransformForward = false for the fixed transform parameters.
   */
  itkSetMacro(TransformForward, bool);

  /** Estimate parameter scales */
  virtual void EstimateScales(ScalesType &scales);

  /** Compute the shift in voxels when deltaParameters is applied onto the
   * current parameters. */
  virtual double ComputeMaximumVoxelShift(ParametersType deltaParameters);

  virtual IndexValueType GetImageDimension() const;

protected:
  OptimizerParameterEstimator();
  ~OptimizerParameterEstimator(){};

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  /** Get the physical coordinates of image corners */
  void SampleWithCornerPoints();

  /** Randomly select some points as samples */
  void SampleImageDomainRandomly();

  /** Set the sample points for computing pixel shifts */
  void SampleImageDomain();

  TransformBase * GetTransform();
  virtual void ComputeJacobianWithRespectToParameters(const VirtualPointType  & p, JacobianType & jacobian) const;

  template< class TContinuousIndexType > void TransformPointToContinuousIndex(
                              const VirtualPointType &point,
                              TContinuousIndexType &mappedIndex);

  void EstimateScalesFromMaximumShift(ScalesType &parameterScales);
  void EstimateScalesFromJacobian(ScalesType &parameterScales);

  /** The templated version of EstimateScalesFromMaximumShift.
   *  The template argument may be either MovingTransformType
   *  or FixedTransformType.
   */
  template <class TTransform> double ComputeTemplatedMaximumVoxelShift(
                              ParametersType deltaParameters);

private:
  OptimizerParameterEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // the metric object
  MetricPointer                 m_Metric;

  // the transform objects
  FixedTransformPointer         m_FixedTransform;
  MovingTransformPointer        m_MovingTransform;

  // the images
  FixedImagePointer             m_FixedImage;
  MovingImagePointer            m_MovingImage;

  // the virtual image for symmetric registration
  VirtualImagePointer           m_VirtualImage;

  // the image samples in the virtual image domain
  std::vector<VirtualPointType> m_ImageSamples;

  /** m_TransformForward specifies which transform scales to be estimated.
   * m_TransformForward = true for the moving transform parameters.
   * m_TransformForward = false for the fixed transform parameters.
   */
  bool m_TransformForward;

  ScaleStrategyType m_ScaleStrategy;

}; //class OptimizerParameterEstimator


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOptimizerParameterEstimator.hxx"
#endif

#endif /* __itkOptimizerParameterEstimator_h */
