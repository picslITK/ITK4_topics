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
#ifndef __itkQuasiNewtonObjectOptimizer_h
#define __itkQuasiNewtonObjectOptimizer_h

#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "itkIntTypes.h"
#include "itkGradientDescentObjectOptimizer.h"
#include "itkOptimizerParameterEstimatorBase.h"
#include <string>

namespace itk
{
/** \class QuasiNewtonObjectOptimizer
 * \brief Implement a Quasi-Newton optimizer with BFGS Hessian estimation.
 *
 * Second order approximation of the cost function is usually more efficient
 * since it estimates the descent or ascent direction more precisely. However,
 * computation of Hessian is usually expensive or unavailable. Alternatively
 * Quasi-Newton methods can estimate a Hessian from the gradients of previous
 * steps. A commonly used Quasi-Newton method is BFGS.
 *
 * Here BFGS method is used to calcute a Quasi-Newton step. The line search
 * is implemented using Strong Wolfe's Conditions.
 *
 * If a helper object of OptimizerParameterEstimatorBase is set, it is used
 * to estimate the maximum step size in line search.
 *
 * The line search algorithm is from "Introduction to Nonlinear Optimization"
 * by Paul J Atzberger, on
 * Page 7 on http://www.math.ucsb.edu/~atzberg/finance/nonlinearOpt.pdf
 *
 * \ingroup ITKHighDimensionalOptimizers
 */
class ITK_EXPORT QuasiNewtonObjectOptimizer:
  public GradientDescentObjectOptimizer
{
public:
  /** Standard class typedefs. */
  typedef QuasiNewtonObjectOptimizer        Self;
  typedef GradientDescentObjectOptimizer    Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(QuasiNewtonObjectOptimizer, AutomaticGradientDescentOptimizer);

  /** Type for Hessian matrix in the Quasi-Newton method */
  typedef itk::Array2D<double>                      HessianType;

  /** Pointer of OptimizerParameterEstimatorBase. */
  typedef OptimizerParameterEstimatorBase::Pointer  OptimizerParameterEstimatorBasePointer;

  /** Connect the OptimizerParameterEstimator .  */
  itkSetObjectMacro(OptimizerParameterEstimator, OptimizerParameterEstimatorBase);
  itkGetObjectMacro(OptimizerParameterEstimator, OptimizerParameterEstimatorBase);

  /** Set the flag for line search */
  itkSetMacro(LineSearchEnabled, bool);

  /** Start and run the optimization */
  virtual void StartOptimization();

  /** Resume the optimization. Can be called after StopOptimization to
   * resume. The bulk of the optimization work loop is here. */
  virtual void ResumeOptimization();

  /** Advance one step following the Quasi-Newton direction. */
  void AdvanceOneStep(void);

  void AdvanceWithBacktrackingLineSearch(ParametersType direction, double maxStepSize);
  double AdvanceWithStrongWolfeLineSearch(ParametersType direction, double maxStepSize);
  double LineSearchZoom(ParametersType initPosition, double f0, double g0, ParametersType direction, double tlow, double thigh);

protected:

  /** The helper object to estimate the learning rate and scales */
  OptimizerParameterEstimatorBasePointer  m_OptimizerParameterEstimator;

  /** m_MaximumVoxelShift is used  to estimated the largest step size */
  double                                  m_MaximumVoxelShift;

  /** m_MinimumGradientNorm is used to stop optimization at a very small gradient */
  double                                  m_MinimumGradientNorm;

  /** m_MinimumGradientNorm is used to stop optimization at a very small value change */
  double                                  m_MinimumValueChange;

  /** Switch for doing line search */
  bool            m_LineSearchEnabled;

  /** A flag to show if GetValueAndDerivate() is evaluated after applying a change.
   *  This happens when line search is done. */
  bool            m_ValueAndDerivateEvaluated;

  /** The information about the previous step */
  double          m_PreviousValue;
  ParametersType  m_PreviousPosition;
  DerivativeType  m_PreviousGradient;

  /** The Quasi-Newton step */
  ParametersType  m_NewtonStep;

  /** The Hessian and its inverse estimated by a Quasi-Newton method */
  HessianType     m_Hessian;
  HessianType     m_HessianInverse;

  //may be used in the future for Levenberg Marquardt
  //double          m_Lambda;

  /** Do line search on the direction of the Newton step */
  //void LineSearch();

  /** Estimate the Newton step that minimizes the local 2nd
   * order approximation */
  void EstimateNewtonStep();

  /** Estimate the Hessian with BFGS method */
  void EstimateHessian();

  /** Reset Hessian to Identity and use the gradient step*/
  void ResetNewtonStep();

  /** Function to compute the learning rate. */
  virtual double EstimateLearningRate(ParametersType step);

  /** This method is used to watch the turning angles between steps. */
  void DebugStepSizeAndAngles(const char *debugLabel,
                              ParametersType lastStep,
                              ParametersType thisStep) const;

  QuasiNewtonObjectOptimizer();
  virtual ~QuasiNewtonObjectOptimizer()
    {
    }
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** helper object to estimate scales and compute maximum voxel shift. */
  //OptimizerHelper::Pointer    m_OptimizerHelper;

private:
  QuasiNewtonObjectOptimizer(const Self &);     //purposely not implemented
  void operator=(const Self &);           //purposely not implemented

};
} // end namespace itk

#endif
