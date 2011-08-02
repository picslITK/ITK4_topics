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
#include "itkOptimizerParameterEstimatorBase.h"
#include "itkImageRandomConstIteratorWithIndex.h"

#include <iostream>

namespace itk
{

/** \class OptimizerParameterEstimator
 *  \brief Implements an optimizer helper class for estimating scales of
 * Transform parameters and computing the maximum voxel shift from a
 * specific transform.
 *
 * Its input include the fixed and moving images and specific transform
 * parameters. These information is usually unavailable from a general
 * optimizer that is not limited to image registration.
 *
 * Currently we implemented two strategies to estimate scales. In the first
 * strategy ScalesFromShift, the scale of a parameter is estimated from the
 * maximum voxel shift yielded from a unit variation of this parameter. The
 * maximization is computed from checking the corners of the image domain.
 * Using corners is enough for affine transformations for this purpose.
 *
 * In the second strategy ScalesFromJacobian, the scale of a parameter is
 * estimated from the averaged squared norm of the Jacobian w.r.t this
 * parameter. The averaging is done over a uniformly random sampling of the
 * image domain.
 *
 * \ingroup ITK-RegistrationCommon
 */
template < class TMetric, class TFixedTransform, class TMovingTransform >
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

  /** Type of the transform to initialize */
  typedef TFixedTransform                           FixedTransformType;
  typedef typename FixedTransformType::Pointer      FixedTransformPointer;

  typedef TMovingTransform                          MovingTransformType;
  typedef typename MovingTransformType::Pointer     MovingTransformPointer;

  typedef TMetric                                   MetricType;
  typedef typename MetricType::Pointer              MetricPointer;

  /** Image Types to use in the initialization of the transform */
  typedef typename TMetric::FixedImageType              FixedImageType;
  typedef typename TMetric::MovingImageType             MovingImageType;
  typedef typename TMetric::VirtualImageType            VirtualImageType;

  typedef typename FixedImageType::ConstPointer   FixedImagePointer;
  typedef typename MovingImageType::ConstPointer  MovingImagePointer;
  typedef typename VirtualImageType::ConstPointer VirtualImagePointer;

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

  itkStaticConstMacro(ImageDimension, unsigned int,
                      FixedImageType::ImageDimension);

  typedef itk::ImageRegion< itkGetStaticConstMacro(ImageDimension) >
                                                              ImageRegionType;
  typedef itk::Size< itkGetStaticConstMacro( ImageDimension ) >   SizeType;

  /** Standard coordinate point type for this class   */
  typedef typename FixedImageType::PointType                      PointType;

  /** Index and Point typedef support. */
  typedef itk::Index< itkGetStaticConstMacro( ImageDimension ) >  IndexType;

  /** Set the metric used in the registration process */
  //itkSetConstObjectMacro(Metric, MetricType);
  virtual void SetMetric(MetricType *metric)
    {
    if (this->m_Metric != metric)
      {
      this->m_Metric = metric;
      this->m_FixedImage   = metric->GetFixedImage();
      this->m_MovingImage  = metric->GetMovingImage();
      this->m_VirtualImage = metric->GetVirtualDomainImage();

      this->Modified();
      }
    }
  /** Get the metric. */
  itkGetConstObjectMacro(Metric, MetricType);

  /** Set the fixed image used in the registration process */
  itkSetConstObjectMacro(FixedImage,  FixedImageType);
  /** Get the fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Set the moving image used in the registration process */
  itkSetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Set the virtual image used in the registration process */
  itkSetConstObjectMacro(VirtualImage, VirtualImageType);
  /** Get the virtual Image. */
  itkGetConstObjectMacro(VirtualImage, VirtualImageType);

  /** Set the transform */
  itkSetObjectMacro(FixedTransform,  FixedTransformType);
  /** Get the transform */
  itkGetObjectMacro(MovingTransform, MovingTransformType);

  /** Set the order of L-norm */
  itkSetMacro(LNorm, int);

  /** Set the flag for forward direction:
   * m_TransformForward = true when the transform mapps from FixedImage
   * domain to MovingImage domain,
   * m_TransformForward = false when the transform mapps from MovingImage
   * domain to FixedImage domain.
   */
  itkSetMacro(TransformForward, bool);

  /** Estimate parameter scales */
  virtual void EstimateScales(ParametersType parameters, ScalesType &scales);

  /** Compute the shift in voxels when deltaParameters is applied onto the
   * current parameters. */
  virtual double ComputeMaximumVoxelShift(ParametersType parameters,
                            ParametersType deltaParameters);

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
  double ComputeLNorm(Point<double, ImageDimension> point);

  /** Set the sample points for computing pixel shifts */
  void SampleImageDomain();

  void EstimateScalesFromMaximumShift(ParametersType parameters,
                                      ScalesType &parameterScales);
  template <class TTransform> void EstimateScalesFromJacobian(ParametersType parameters,
                                  ScalesType &parameterScales);
  template <class TTransform> double ComputeTemplatedMaximumVoxelShift(
                            ParametersType parameters,
                            ParametersType deltaParameters);

private:
  OptimizerParameterEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MetricPointer             m_Metric;

  FixedTransformPointer     m_FixedTransform;
  MovingTransformPointer    m_MovingTransform;

  FixedImagePointer   m_FixedImage;
  MovingImagePointer  m_MovingImage;
  VirtualImagePointer m_VirtualImage;

  std::vector<PointType> m_ImageSamples;

  /** Specify how to calculate the distance between two points */
  int m_LNorm;

  /** Specify the transformation direction. Set to true when the transform
   * mapps from FixedImage domain to MovingImage domain*/
  bool m_TransformForward;

  ScaleStrategyType m_ScaleStrategy;

}; //class OptimizerParameterEstimator


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOptimizerParameterEstimator.hxx"
#endif

#endif /* __itkOptimizerParameterEstimator_h */
