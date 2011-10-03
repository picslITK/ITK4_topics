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
#ifndef __itkRegistrationParameterScalesFromShift_h
#define __itkRegistrationParameterScalesFromShift_h

#include "itkRegistrationParameterScalesEstimator.h"

namespace itk
{

/** \class RegistrationParameterScalesFromShift
 *  \brief Implements a registration helper class for estimating scales of
 * transform parameters from the maximum voxel shift caused by a parameter
 * change.
 *
 * Its input includes the fixed/moving images and transform objects,
 * which can be obtained from the metric object.
 *
 * The scale of a parameter is estimated from the maximum voxel shift produced
 * from a small variation of this parameter. The maximization is done by checking
 * the pixels of the image domain. By defaults, all pixels will be checked.
 * For affine transforms, we may check only the corners by calling
 * SetSamplingStrategy(CornerSampling).
 *
 * \ingroup ITKHighDimensionalOptimizers
 */
template < class TMetric >
class ITK_EXPORT RegistrationParameterScalesFromShift :
  public RegistrationParameterScalesEstimator< TMetric >
{
public:
  /** Standard class typedefs. */
  typedef RegistrationParameterScalesFromShift            Self;
  typedef RegistrationParameterScalesEstimator< TMetric > Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( RegistrationParameterScalesFromShift, RegistrationParameterScalesEstimator );

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

  /** Estimate parameter scales */
  virtual void EstimateScales(ScalesType &scales);

protected:
  RegistrationParameterScalesFromShift();
  ~RegistrationParameterScalesFromShift(){};

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  /** Compute the shift in voxels when deltaParameters is applied onto the
   * current parameters. */
  virtual FloatType ComputeMaximumVoxelShift(const ParametersType &deltaParameters);

  /** The templated version of EstimateScalesFromMaximumShift.
   *  The template argument may be either MovingTransformType
   *  or FixedTransformType.
   */
  template <class TTransform> FloatType ComputeMaximumVoxelShiftTemplated(
                              const ParametersType &deltaParameters);

private:
  RegistrationParameterScalesFromShift(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //A small variation of parameters
  typename ParametersType::ValueType   m_SmallParameterVariation;

}; //class RegistrationParameterScalesFromShift


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRegistrationParameterScalesFromShift.hxx"
#endif

#endif /* __itkRegistrationParameterScalesFromShift_h */
