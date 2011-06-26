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
#ifndef __itkQuasiNewtonOptimizer_h
#define __itkQuasiNewtonOptimizer_h

#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "itkIntTypes.h"
#include "itkAutomaticGradientDescentOptimizer.h"
#include "itkOptimizerHelper.h"
#include <string>
namespace itk
{
/** \class QuasiNewtonOptimizer
 * \brief Implement a Quasi-Newton optimizer with BFGS Hessian estimation.
 *
 * Second order approximation of the cost function is usually more efficient
 * since it estimates the descent or ascent step more precisely. However, the
 * computation of Hessian is usually expensive or unavailable. Alternatively
 * Quasi-Newton methods can estimate a Hessian from the gradients of previous
 * steps. A commonly used Quasi-Newton method is BFGS.
 *
 * It uses BFGS method to calcute a Quasi-Newton step. The backtracking line
 * search is implemented but disabled by default. Instead, the voxel shift
 * is used in restricting step sizes.
 *
 * Currently the Hessian is estimated in the scaled parameter space. Since 2nd
 * order step is robust against scales, we will change it to the original
 * parameter space later. In addition, Sherman-Morrison formula may be used
 * for faster computation of the inverse Hessian.
 *
 * \ingroup Numerics Optimizers
 * \ingroup ITK-Optimizers
 */
class ITK_EXPORT QuasiNewtonOptimizer:
  public AutomaticGradientDescentOptimizer
{
public:
  /** Standard class typedefs. */
  typedef QuasiNewtonOptimizer              Self;
  typedef AutomaticGradientDescentOptimizer Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(QuasiNewtonOptimizer, AutomaticGradientDescentOptimizer);

  /** Type for Hessian matrix in the Quasi-Newton method */
  typedef Array2D<double>            HessianType;

  /** Advance one step following the gradient direction. */
  void AdvanceOneStep(void);

protected:

  /** The gradient in the previous step */
  DerivativeType  m_PreviousGradient;
  ParametersType  m_PreviousPosition;

  /** The hessian estimated by a Quasi-Newton method
   */
  HessianType     m_Hessian;
  HessianType     m_HessianInverse;
  DerivativeType  m_NewtonStep;
  DerivativeType  m_ScaledNewtonStep;
  DerivativeType  m_PreviousNewtonStep;

  bool            m_LineSearchEnabled;

  /** Do line search on the direction of the Newton step */
  void LineSearch();

  /** Estimate the Newton step that minimizes the local 2nd
   * order approximation */
  void EstimateNewtonStep();

  /** Estimate the Hessian with BFGS method */
  void EstimateHessian(ParametersType x1, ParametersType x2,
                       DerivativeType g1, DerivativeType g2);

  /** Translate the parameters into the scaled space */
  void ScalePosition(ParametersType p1, ParametersType &p2);

  /** Translate the parameters into the original space */
  void ScaleBackPosition(ParametersType p1, ParametersType &p2);

  /** Translate the derivative into the scaled space */
  void ScaleDerivative(DerivativeType g1, DerivativeType &p2);

  /** Translate the derivative into the original space */
  void ScaleBackDerivative(DerivativeType g1, DerivativeType &p2);

  QuasiNewtonOptimizer();
  virtual ~QuasiNewtonOptimizer() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  QuasiNewtonOptimizer(const Self &);     //purposely not implemented
  void operator=(const Self &);           //purposely not implemented

};
} // end namespace itk

#endif
