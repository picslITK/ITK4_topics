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
 * We've implemented two strategies to estimate scales. With the first
 * strategy ScalesFromShift, the scale of a parameter is estimated from the
 * maximum voxel shift yielded from a unit variation of this parameter. The
 * maximization is done by checking the corners of the image domain.
 * Using corners is enough for affine transformations for this purpose, but
 * may be changed for other types of transforms.
 *
 * With the second strategy ScalesFromJacobian, the scale of a parameter is
 * estimated from the averaged squared norm of the Jacobian w.r.t this
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

  /** Jacobian type */
  typedef typename Superclass::JacobianType         JacobianType;

  /** Image Types to use in the initialization of the transform */
  typedef typename TMetric::FixedImageType          FixedImageType;
  typedef typename TMetric::MovingImageType         MovingImageType;
  typedef typename TMetric::VirtualImageType        VirtualImageType;

  typedef typename FixedImageType::ConstPointer     FixedImagePointer;
  typedef typename MovingImageType::ConstPointer    MovingImagePointer;
  typedef typename VirtualImageType::ConstPointer   VirtualImagePointer;

  /* Image dimension accessors */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, unsigned int,
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

  /** Set the metric used in the registration process */
  virtual void SetMetric(MetricType *metric)
    {
    if (this->m_Metric != metric)
      {
      this->m_Metric = metric;
      this->m_FixedImage   = metric->GetFixedImage();
      this->m_MovingImage  = metric->GetMovingImage();
      this->m_VirtualImage = metric->GetVirtualDomainImage();

      SetMovingTransform(const_cast<MovingTransformType *>(m_Metric->GetMovingTransform()));
      SetFixedTransform(const_cast<FixedTransformType *>(m_Metric->GetFixedTransform()));

      this->Modified();
      }
    }
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
  /** Set the moving transform */
  itkSetObjectMacro(MovingTransform, MovingTransformType);

  /** Set the order of L-norm for computing point distances */
  itkSetMacro(LNorm, int);

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

  virtual unsigned int GetImageDimension()
    {
    return Self::VirtualImageDimension;
    }

protected:
  OptimizerParameterEstimator();
  ~OptimizerParameterEstimator(){};

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  /** Get the physical coordinates of image corners */
  void SampleWithCornerPoints();

  /** Randomly select some points as samples */
  void SampleImageDomainRandomly();

  /** Compute the L-norm of a point */
  template< class TContinuousIndexType > double ComputeLNorm(TContinuousIndexType point);

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

  /** Specify how to calculate the distance between two points */
  int m_LNorm;

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
