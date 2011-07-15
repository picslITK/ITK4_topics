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
#ifndef __itkOptimizerHelper_h
#define __itkOptimizerHelper_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>

namespace itk
{

/** \class OptimizerHelper
 *  \brief OptimizerHelper is a helper class intended to estimate the scales
 * and the maximum voxel shift in optimizers of image registration.
 *
 * The estimation requires information about images and transform that are
 * usually not available in general optimizers that are not limited to
 * image registration. Therefore we create a separate helper class for the
 * specific type of optimization.
 *
 * This is an abstract class to be implemented with different estimation
 * strategies. An example of subclass is itkOptimizerParameterEstimator.
 *
 * \ingroup Numerics Optimizers
 * \ingroup ITK-Optimizers
 */
class ITK_EXPORT OptimizerHelper : public Object
{
public:
  /** Standard class typedefs. */
  typedef OptimizerHelper                       Self;
  typedef Object                                Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  //itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimizerHelper, Object );

  /** Type of scales */
  typedef Array< double > ScalesType;
  /** Type of paramters in optimizer */
  typedef Array< double > ParametersType;
  /** Type of derivative of metric cost function */
  typedef Array< double > DerivativeType;

  /** Estimate parameter scales and learning rate*/
  virtual void EstimateScales(ParametersType parameters,
    ScalesType &scales) = 0;

  /** Compute the maximum voxel shift from the parameter change */
  virtual double ComputeMaximumVoxelShift(ParametersType parameters,
    ParametersType deltaParameters) = 0;

protected:
  OptimizerHelper(){};
  ~OptimizerHelper(){};

  void PrintSelf(std::ostream &os, Indent indent) const
    {
    Superclass::PrintSelf(os,indent);
    }

private:
  OptimizerHelper(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

}; //class OptimizerHelper

}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
//#include "itkOptimizerHelper.txx"
#endif

#endif /* __itkOptimizerHelper_h */
