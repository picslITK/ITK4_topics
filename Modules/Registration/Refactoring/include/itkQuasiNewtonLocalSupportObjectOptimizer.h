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
#ifndef __itkQuasiNewtonLocalSupportObjectOptimizer_h
#define __itkQuasiNewtonLocalSupportObjectOptimizer_h

#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "itkIntTypes.h"
#include "itkGradientDescentObjectOptimizer.h"
//#include "itkImageToImageObjectMetric.h"
#include "itkOptimizerParameterEstimatorBase.h"
#include <string>

namespace itk
{
/** \class QuasiNewtonLocalSupportObjectOptimizer
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
 * For metric classes using a deformation field, we have "local support" to
 * compute the local update of the deformation field. In this case, we compute
 * the local Quasi-Newton step at a specific index of the field.
 *
 * \ingroup ITKRegistrationRefactoring
 */
class ITK_EXPORT QuasiNewtonLocalSupportObjectOptimizer:
  public GradientDescentObjectOptimizer
{
public:
  /** Standard class typedefs. */
  typedef QuasiNewtonLocalSupportObjectOptimizer        Self;
  typedef GradientDescentObjectOptimizer                Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(QuasiNewtonLocalSupportObjectOptimizer, AutomaticGradientDescentOptimizer);

  /** Metric function type */
  //typedef ImageToImageObjectMetric                  ImageMetricType;
  //typedef MetricType::Pointer                       ImageMetricPointer;

  /** Type for Hessian matrix in the Quasi-Newton method */
  typedef itk::Array2D<double>                      HessianType;
  typedef itk::Array2D<double>                      LocalHessianType;

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

  /** Advance one step following the Quasi-Newton direction
   * if the transform has local support. */
  void AdvanceOneLocalStep(void);

protected:

  /** The helper object to estimate the learning rate and scales */
  OptimizerParameterEstimatorBasePointer  m_OptimizerParameterEstimator;
  double                                  m_MaximumVoxelShift;
  double                                  m_MinimumGradientNorm;
  double                                  m_MinimumValueChange;

  /** Switch for doing line search */
  bool            m_LineSearchEnabled;

  /** A flag to show if GetValueAndDerivate() is evaluated after applying a change.
   *  This happens when line search is done. */
  bool            m_ValueAndDerivateEvaluated;

  /** The gradient in the previous step */
  ParametersType  m_PreviousPosition;
  DerivativeType  m_PreviousGradient;

  /** The Quasi-Newton step */
  ParametersType  m_NewtonStep;

  /** The Hessian estimated by a Quasi-Newton method */
  HessianType     m_Hessian;
  HessianType     m_HessianInverse;

  /** The Hessian with local support */
  LocalHessianType  *m_LocalHessian;
  LocalHessianType  *m_LocalHessianInverse;

  /** Do line search on the direction of the Newton step */
  //void LineSearch();

  /** Estimate the Newton step that minimizes the local 2nd
   * order approximation with local support. */
  void EstimateLocalNewtonStep();

  /** Estimate the Hessian with BFGS method with local support. */
  void EstimateLocalHessian();

  /** Function to compute the learning rate. */
  virtual double EstimateLearningRate(ParametersType step);

  /** This method is used to watch the turning angles between steps. */
  void DebugStepSizeAndAngles(const char *debugLabel,
                              ParametersType lastStep,
                              ParametersType thisStep) const;

  QuasiNewtonLocalSupportObjectOptimizer();
  virtual ~QuasiNewtonLocalSupportObjectOptimizer()
    {
    if (m_LocalHessian)
      {
      delete[] m_LocalHessian;
      }
    if (m_LocalHessianInverse)
      {
      delete[] m_LocalHessianInverse;
      }
    }
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** helper object to estimate scales and compute maximum voxel shift. */
  //OptimizerHelper::Pointer    m_OptimizerHelper;

private:
  QuasiNewtonLocalSupportObjectOptimizer(const Self &);     //purposely not implemented
  void operator=(const Self &);           //purposely not implemented

};
} // end namespace itk

#endif
