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
#ifndef __itkParameterScaleEstimator_h
#define __itkParameterScaleEstimator_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>

namespace itk
{

/** \class ParameterScaleEstimator
 *  \brief ParameterScaleEstimator is a helper class intended to
 * initialize the scales of Transform parameters in registration.
 *
 * This class is connected to the fixed image, moving image and transform
 * involved in the registration.
 *
 * Define a transform as Y = T(P,X) where X, Y are coordinates, and P are
 * parameters. Let p be one of the parameters. Introduce an intermediate
 * variable q = s * p where s is the scaling factor. The p column of the
 * transform Jacobian is written as d_Y/d_p.
 *
 * The scales are estimated in a way such that a unit change of each
 * parameter would yield roughly similar voxel shifts. That is, d_Y/d_q
 * is roughly same for all scaled parameters.
 *
 * After applying scale s, d_Y/d_q = d_Y/d_p * (1/s). Therefore we may
 * choose s = d_Y/d_p.

 * And we have the metric derivative d_metric/d_q = d_metric/d_p * (1/s).
 *
 * Let us consider the impact of s to the scales in ITK optimizers. In the
 * ITK class GradientDescentOptimizer, the step of parameter p is
 *   step_p = (d_metric/d_p) / old_opt_scale * learningRate .
 *
 * With q = s * p, the step for parameter q will be
 *   step_q = (d_metric/d_q) / old_opt_scale * learningRate
 *          = (d_metric/d_p) / s / old_opt_scale * learningRate .
 *
 * Since step_p = step_q / s, we have
 *   step_p = (d_metric/d_p) / s / s / old_opt_scale * learningRate .
 *
 * With parameter scale s, the optimizer scale for p is changed to
 *   new_opt_scale = s * s * old_opt_scale .
 * where s = d_Y/d_p.
 *
 * \ingroup ITK-Registration
 */
template < class TFixedImage,
           class TMovingImage,
           class TTransform >
class ITK_EXPORT ParameterScaleEstimator : public Object
{
public:
  /** Standard class typedefs. */
  typedef ParameterScaleEstimator               Self;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ParameterScaleEstimator, Object );

  /** Type of the transform to initialize */
  typedef TTransform                                TransformType;
  typedef typename TransformType::Pointer           TransformPointer;
  typedef typename TransformType::ParametersType    ParametersType;

  /** Scale type.
   * This array defines scale to be applied to parameters before
   * being evaluated in the cost function. This allows to
   * map to a more convenient space. In particular this is
   * used to normalize parameter spaces in which some parameters
   * have a different dynamic range.   */
  typedef Array< double > ScalesType;
  typedef Array< double > DerivativeType;

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
  itkSetConstObjectMacro( FixedImage,  FixedImageType  );

  /** Set the moving image used in the registration process */
  itkSetConstObjectMacro( MovingImage, MovingImageType );

  /** Set the order of L-norm */
  itkSetMacro(LNorm, int);

  /** Set the global scaling factor */
  itkSetMacro(GlobalScalingFactor, double);

  /** Set the flag for forward direction:
   * m_TransformForward = true when the transform mapps from FixedImage
   * domain to MovingImage domain,
   * m_TransformForward = false when the transform mapps from MovingImage
   * domain to FixedImage domain.
   */
  itkSetMacro(TransformForward, bool);

  /** Estimate parameter scales */
  void EstimateScales(TransformPointer transform, ScalesType &parameterScales);

  /** Estimate parameter scales and learningRate*/
  void EstimateScalesAndLearningRate(TransformPointer transform,
                                DerivativeType   gradient,
                                ScalesType       &parameterScales,
                                double           &learningRate);

  /** Set the sample points for computing pixel shifts */
  void SampleImageDomain();

protected:
  ParameterScaleEstimator();
  ~ParameterScaleEstimator(){};

  void PrintSelf(std::ostream &os, Indent indent) const;

  /** Get the physical coordinates of image corners */
  template <class ImageType>
  void GetImageCornerPoints(const ImageType *image,
                            std::vector<PointType> &samples);

  /** Compute the maximum shift when deltaParameters is applied onto the
   * current transform. */
  void ComputeMaximumShift(typename TransformType::Pointer transform,
                           ParametersType deltaParameters, double &maxShift);

  /** Compute the L-norm of a point */
  double ComputeLNorm(Point<double, ImageDimension> point);

private:
  ParameterScaleEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  FixedImagePointer   m_FixedImage;
  MovingImagePointer  m_MovingImage;

  std::vector<PointType> m_FixedImageSamples;
  std::vector<PointType> m_MovingImageSamples;

  /** Specify how to calculate the distance between two points */
  int m_LNorm;

  /** Change all scales by this denominator */
  double m_GlobalScalingFactor;

  /** Specify the transformation direction. Set to true when the transform
   * mapps from FixedImage domain to MovingImage domain*/
  bool m_TransformForward;
}; //class ParameterScaleEstimator


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParameterScaleEstimator.txx"
#endif

#endif /* __itkParameterScaleEstimator_h */
