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
#ifndef __itkOptimizerParameterEstimatorBase_h
#define __itkOptimizerParameterEstimatorBase_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkArray.h"
#include "itkArray2D.h"

#include <iostream>

namespace itk
{

/** \class OptimizerParameterEstimatorBase
 *  \brief OptimizerParameterEstimatorBase is the base class to offer some
 * empty methods of estimating the scales and the maximum voxel shift for
 * optimizers of image registration.
 *
 * The estimation requires information about images and transform objects
 * that are usually not available in general optimizers. Therefore we create
 * this helper class for optimizers about image registration.
 *
 * A default implementation of OptimizerParameterEstimatorBase is in
 * itkOptimizerParameterEstimator.
 *
 * \ingroup ITKHighDimensionalOptimizers
 */
class ITK_EXPORT OptimizerParameterEstimatorBase : public Object
{
public:
  /** Standard class typedefs. */
  typedef OptimizerParameterEstimatorBase       Self;
  typedef Object                                Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimizerParameterEstimatorBase, Object );

  /** Type of scales */
  typedef Array< double >     ScalesType;
  /** Type of paramters in optimizer */
  typedef Array< double >     ParametersType;
  /** Type of derivative of metric cost function */
  typedef Array< double >     DerivativeType;
  /** Type of Jacobian of transform */
  typedef Array2D< double >   JacobianType;

  /** Estimate parameter scales and learning rate*/
  virtual void EstimateScales(ScalesType &scales) = 0;

  /** Compute the maximum voxel shift from the parameter change */
  virtual double ComputeMaximumVoxelShift(ParametersType deltaParameters) = 0;

  virtual unsigned int GetImageDimension() = 0;

protected:
  OptimizerParameterEstimatorBase(){};
  ~OptimizerParameterEstimatorBase(){};

  void PrintSelf(std::ostream &os, Indent indent) const
    {
    Superclass::PrintSelf(os,indent);
    }

private:
  OptimizerParameterEstimatorBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

}; //class OptimizerParameterEstimatorBase

}  // namespace itk

#endif /* __itkOptimizerParameterEstimatorBase_h */
