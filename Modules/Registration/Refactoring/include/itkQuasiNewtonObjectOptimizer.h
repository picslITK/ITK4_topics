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
#include "itkImageToImageObjectMetric.h"
#include <string>

namespace itk
{
/** \class QuasiNewtonObjectOptimizer
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

  /** Metric function type */
  //typedef ImageToImageObjectMetric                  ImageMetricType;
  //typedef MetricType::Pointer                       ImageMetricPointer;

  /** Type for Hessian matrix in the Quasi-Newton method */
  typedef Array2D<double>            HessianType;

  /** Methods to configure the cost function. */
  itkGetConstReferenceMacro(Maximize, bool);
  itkSetMacro(Maximize, bool);
  itkBooleanMacro(Maximize);
  bool GetMinimize() const
  { return !m_Maximize; }
  void SetMinimize(bool v)
  { this->SetMaximize(!v); }
  void MinimizeOn()
  { this->MaximizeOff(); }
  void MinimizeOff()
  { this->MaximizeOn(); }

  /** Start and run the optimization */
  virtual void StartOptimization();

  /** Advance one step following the Quasi-Newton direction. */
  void AdvanceOneStep(void);

  /** Advance one step following the Quasi-Newton direction
   * if the transform has local support. */
  void AdvanceOneLocalStep(void);

protected:

  /** The gradient in the previous step */
  ParametersType  m_PreviousPosition;
  DerivativeType  m_PreviousGradient;

  // Reference to current position
  // Use reference to save memory copy
  //ParametersType& m_CurrentPositionRef;

  /** The hessian estimated by a Quasi-Newton method
   */
  HessianType     m_Hessian;
  HessianType     m_HessianInverse;

  ParametersType  m_LocalHessian;
  ParametersType  m_NewtonStep;
  ParametersType  m_ScaledNewtonStep;

  //bool            m_LineSearchEnabled;

  /** Do line search on the direction of the Newton step */
  //void LineSearch();

  /** Estimate the Newton step that minimizes the local 2nd
   * order approximation */
  void EstimateNewtonStep();

  /** Estimate the Hessian with BFGS method */
  void EstimateHessian();

  /** Estimate the Newton step that minimizes the local 2nd
   * order approximation with local support. */
  void EstimateLocalNewtonStep();

  /** Estimate the Hessian with BFGS method with local support. */
  void EstimateLocalHessian();

  /** Translate the parameters into the scaled space */
  void ScalePosition(ParametersType p1, ParametersType &p2);

  /** Translate the parameters into the original space */
  void ScaleBackPosition(ParametersType p1, ParametersType &p2);

  /** Translate the derivative into the scaled space */
  void ScaleDerivative(DerivativeType g1, DerivativeType &p2);

  /** Translate the derivative into the original space */
  void ScaleBackDerivative(DerivativeType g1, DerivativeType &p2);

  double                        m_MaximumVoxelShift;
  double                        m_MinimumVoxelShift;

  StopConditionType             m_StopCondition;
  std::ostringstream            m_StopConditionDescription;

  /** Function to compute the learning rate. */
  virtual double EstimateLearningRate(ParametersType step);

  /** This method is used to watch the turning angles between steps. */
  void DebugStepSizeAndAngles(const char *debugLabel,
                              ParametersType lastStep,
                              ParametersType thisStep) const;

  QuasiNewtonObjectOptimizer();
  virtual ~QuasiNewtonObjectOptimizer() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  bool    m_Maximize;
  double  m_LearningRate;
  bool    m_HasLocalSupport;

private:
  QuasiNewtonObjectOptimizer(const Self &);     //purposely not implemented
  void operator=(const Self &);           //purposely not implemented

};
} // end namespace itk

#endif
