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
#include "itkOptimizerHelper.h"

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
 * \ingroup ITK-Registration
 */
template < class TFixedImage,
           class TMovingImage,
           class TTransform >
class ITK_EXPORT OptimizerParameterEstimator : public OptimizerHelper
{
public:
  /** Standard class typedefs. */
  typedef OptimizerParameterEstimator           Self;
  typedef OptimizerHelper                       Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimizerParameterEstimator, OptimizerHelper );

  /** Type of the transform to initialize */
  typedef TTransform                                TransformType;
  typedef typename TransformType::Pointer           TransformPointer;
  typedef typename TransformType::ParametersType    ParametersType;

  /** The stratigies to decide scales */
  typedef enum { ScalesFromShift, ScalesFromJacobian } ScaleStrategyType;
  /** Set the learning rate strategy */
  itkSetMacro(ScaleStrategy, ScaleStrategyType);

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TFixedImage::ImageDimension);

  typedef itk::ImageRegion< itkGetStaticConstMacro(ImageDimension) >
                                                              ImageRegionType;
  typedef itk::Size< itkGetStaticConstMacro( ImageDimension ) >   SizeType;

  /** Standard coordinate point type for this class   */
  typedef typename TFixedImage::PointType                         PointType;

  /** Index and Point typedef support. */
  typedef itk::Index< itkGetStaticConstMacro( ImageDimension ) >  IndexType;

  /** Image Types to use in the initialization of the transform */
  typedef   TFixedImage              FixedImageType;
  typedef   TMovingImage             MovingImageType;

  typedef   typename FixedImageType::ConstPointer   FixedImagePointer;
  typedef   typename MovingImageType::ConstPointer  MovingImagePointer;

  /** Set the fixed image used in the registration process */
  itkSetConstObjectMacro(FixedImage,  FixedImageType);
  /** Get the Fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Set the moving image used in the registration process */
  itkSetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the Moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Set the transform */
  itkSetObjectMacro(Transform, TransformType);
  /** Get the transform */
  itkGetObjectMacro(Transform, TransformType);

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
  void EstimateScales(ParametersType parameters, ScalesType &scales);

  /** Compute the shift in voxels when deltaParameters is applied onto the
   * current parameters. */
  double ComputeMaximumVoxelShift(ParametersType parameters,
                           ParametersType deltaParameters);


protected:
  OptimizerParameterEstimator();
  ~OptimizerParameterEstimator(){};

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  /** Get the physical coordinates of image corners */
  template <class ImageType> void SampleWithCornerPoints();

  /** Randomly select some points as samples */
  template <class ImageType> void SampleImageDomainRandomly();

  /** Compute the L-norm of a point */
  double ComputeLNorm(Point<double, ImageDimension> point);

  /** Set the sample points for computing pixel shifts */
  void SampleImageDomain();

  void EstimateScalesFromMaximumShift(ParametersType parameters,
                                      ScalesType &parameterScales);
  void EstimateScalesFromJacobian(ParametersType parameters,
                                  ScalesType &parameterScales);

private:
  OptimizerParameterEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  FixedImagePointer   m_FixedImage;
  MovingImagePointer  m_MovingImage;
  TransformPointer    m_Transform;

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
#include "itkOptimizerParameterEstimator.txx"
#endif

#endif /* __itkOptimizerParameterEstimator_h */
