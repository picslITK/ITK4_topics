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
#ifndef _itkQuasiNewtonObjectOptimizer_txx
#define _itkQuasiNewtonObjectOptimizer_txx

#include "itkQuasiNewtonObjectOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkMacro.h"

namespace itk
{
/**
 * Constructor
 */
QuasiNewtonObjectOptimizer
::QuasiNewtonObjectOptimizer()
{
  itkDebugMacro("Constructor");

  m_MaximumVoxelShift = 1.0;

  // m_MinimumVoxelShift will be set to a value in the first iteration
  // according the first step
  m_MinimumVoxelShift = 0;

}

void
QuasiNewtonObjectOptimizer
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * Start the optimization
 */
void
QuasiNewtonObjectOptimizer
::StartOptimization(void)
{
  itkDebugMacro("StartOptimization");

  if (!this->m_Metric->HasLocalSupport())
    {
    // initialize scales
    ScalesType scales(this->m_Metric->GetNumberOfParameters());
    this->m_Metric->EstimateScales(true, scales);
    //m_OptimizerHelper->EstimateScales(this->m_Metric->GetParameters(), scales);
    this->SetScales(scales);
    std::cout << " Estimated scales = " << scales << std::endl;
    }

  this->Superclass::StartOptimization();

}

/**
 * Advance one Step: if the Quasi-Newton direction is consistent
 * with the gradient direction, follow the Quasi-Newton direction,
 * otherwise, follow the gradient direction.
 */
void
QuasiNewtonObjectOptimizer
::AdvanceOneStep(void)
{
  itkDebugMacro("AdvanceOneStep");

  if ( this->m_Metric->HasLocalSupport() )
    {
    AdvanceOneLocalStep();
    return;
    }

  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  ScalesType scales = this->GetScales();

  double learningRate;
  bool   newtonStepException = false;

  this->m_CurrentPosition = this->m_Metric->GetParameters();

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousPosition = this->m_CurrentPosition;
    m_PreviousGradient = this->m_Gradient;
    }

  try
    {
    this->EstimateNewtonStep();
    }
  catch ( ExceptionObject & )
    {
    //This may happen with a singular hessian matrix
    newtonStepException = true;
    }

  /** Save for the next iteration */
  m_PreviousPosition = this->GetCurrentPosition();
  m_PreviousGradient = this->GetGradient();

  DerivativeType gradientStep = m_Gradient;
  for (unsigned int i=0; i<spaceDimension; i++)
    {
    gradientStep[i] = gradientStep[i] / scales[i];
    }

  /** If a Newton step is on the opposite direction of a gradient step, we'd
   * better use the gradient step. This happens when the second order
   * approximation produces a convex instead of an expected concave, or
   * vice versa.
   */
  if ( newtonStepException || inner_product(gradientStep, m_NewtonStep) <= 0 )
    {
    learningRate = this->EstimateLearningRate(gradientStep);
    this->SetLearningRate( learningRate );

    if ( learningRate == 0)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Learning rate is zero after "
                                 << this->GetCurrentIteration()
                                 << " iterations. This may be due to that"
                                 << " the new step yields zero voxel shift.";
      this->StopOptimization();
      return;
      }
    this->GradientDescentObjectOptimizer::AdvanceOneStep();
    }
  else
    {
    // Now a Newton step is on the consistent direction of a gradient step
    learningRate = this->EstimateLearningRate(m_NewtonStep);
    learningRate = vnl_math_min(learningRate, 1.0);
    this->SetLearningRate( learningRate );

    if ( learningRate == 0)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Learning rate is zero after "
                                 << this->GetCurrentIteration()
                                 << " iterations. This may be due to that"
                                 << " the new step yields little voxel shift.";
      this->StopOptimization();
      return;
      }

    DerivativeType step(spaceDimension);
    for ( unsigned int j = 0; j < spaceDimension; j++ )
      {
      step[j] = this->m_LearningRate * this->m_NewtonStep[j];
      }
    this->m_Metric->UpdateTransformParameters( step );

    this->InvokeEvent( IterationEvent() );
    }
}

/**
 * Advance one Step: if the Quasi-Newton direction is consistent
 * with the gradient direction, follow the Quasi-Newton direction,
 * otherwise, follow the gradient direction.
 */
void
QuasiNewtonObjectOptimizer
::AdvanceOneLocalStep(void)
{
  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  ScalesType scales = this->GetScales();

  // Use reference to save memory copy
  const ParametersType & currentPositionRef = this->m_Metric->GetParameters();

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousPosition = currentPositionRef;
    m_PreviousGradient = this->m_Gradient;
    }

  try
    {
    this->EstimateLocalNewtonStep();
    }
  catch ( ExceptionObject & excp )
    {
    m_StopCondition = QuasiNewtonStepError;
    m_StopConditionDescription << "QuasiNewton step error after "
                               << this->GetCurrentIteration()
                               << " iterations. "
                               << excp.GetDescription();
    this->StopOptimization();
    return;
    }

  /** Save for the next iteration */
  m_PreviousPosition = currentPositionRef;
  m_PreviousGradient = this->GetGradient();

  double maxGradient = 0, maxNewtonStep = 0;
  for ( unsigned int p=0; p<spaceDimension; p++ )
    {
    if (maxGradient < vcl_abs(m_Gradient[p]))
      {
      maxGradient = vcl_abs(m_Gradient[p]);
      }
    if (maxNewtonStep < vcl_abs(m_NewtonStep[p]))
      {
      maxNewtonStep = vcl_abs(m_NewtonStep[p]);
      }
    }

  double scaleGradient = maxGradient / m_MaximumVoxelShift;
  double scaleNewtonStep = maxNewtonStep / m_MaximumVoxelShift;

  for ( unsigned int p=0; p<spaceDimension; p++ )
    {
    /** If a Newton step is on the opposite direction of a gradient step, we'd
     * better use the gradient step. This happens when the second order
     * approximation produces a convex instead of an expected concave, or
     * vice versa.
     */
    if ( (m_Gradient[p] * m_NewtonStep[p]) <= 0 )
      {
      m_NewtonStep[p] = m_Gradient[p] / scaleGradient;
      }
    else
      {
      m_NewtonStep[p] = m_NewtonStep[p] / scaleNewtonStep;
      }

    } //end of for

  this->m_Metric->UpdateTransformParameters( this->m_NewtonStep );

  this->InvokeEvent( IterationEvent() );
}

/** Estimate Hessian step */
void QuasiNewtonObjectOptimizer
::EstimateNewtonStep()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  // Estimate Hessian
  EstimateHessian();

  // Compute the Newton step
  ParametersType sdx(numPara);
  sdx = m_HessianInverse * m_Gradient;

  this->m_NewtonStep = -1.0 * sdx;

}

/** Estimate Hessian matrix */
void QuasiNewtonObjectOptimizer
::EstimateHessian()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  // Initialize Hessian to identity matrix
  if ( this->GetCurrentIteration() == 0 )
    {
    m_Hessian.SetSize(numPara, numPara);
    m_Hessian.Fill(0.0f);
    m_HessianInverse.SetSize(numPara, numPara);
    m_HessianInverse.Fill(0.0f);

    for (unsigned int i=0; i<numPara; i++)
      {
      m_Hessian[i][i] = 1.0; //identity matrix
      m_HessianInverse[i][i] = 1.0; //identity matrix
      }
    return;
    }

  ParametersType dx(numPara);  //delta of position x: x_k+1 - x_k
  ParametersType dg(numPara);  //delta of gradient: g_k+1 - g_k
  ParametersType edg(numPara); //estimated delta of gradient: hessian_k * dx

  dx = this->m_CurrentPosition - this->m_PreviousPosition;
  dg = this->m_Gradient - this->m_PreviousGradient;
  edg = m_Hessian * dx;

  double dot_dg_dx = inner_product(dg, dx);
  double dot_edg_dx = inner_product(edg, dx);

  if (dot_dg_dx ==0 || dot_edg_dx == 0)
    {
    itkExceptionMacro(<< "Division by zero in Quasi-Newton step. ");
    return;
    }

  vnl_matrix<double> plus  = outer_product(dg, dg) / dot_dg_dx;
  vnl_matrix<double> minus = outer_product(edg, edg) / dot_edg_dx;
  vnl_matrix<double> newHessian = m_Hessian + plus - minus;

  m_Hessian         = newHessian;
  m_HessianInverse  = vnl_matrix_inverse<double>(newHessian);

}

/** Estimate Hessian step */
void QuasiNewtonObjectOptimizer
::EstimateLocalNewtonStep()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  // Estimate Hessian
  EstimateLocalHessian();

  if ( this->GetCurrentIteration() == 0 )
    {
    m_NewtonStep.SetSize(numPara);
    }

  // Compute the Newton step
  for (unsigned int i=0; i<numPara; i++)
    {
    m_NewtonStep[i] = - m_Gradient[i] / m_LocalHessian[i];
    }

}

/** Estimate Hessian matrix */
void QuasiNewtonObjectOptimizer
::EstimateLocalHessian()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();
  // Use reference to save memory copy
  const ParametersType & currentPositionRef = this->m_Metric->GetParameters();

  // Initialize Hessian to identity matrix
  if ( this->GetCurrentIteration() == 0 )
    {
    m_LocalHessian.SetSize(numPara);
    m_LocalHessian.Fill(1.0f);
    //m_LocalHessianInverse.SetSize(numPara);
    //m_LocalHessianInverse.Fill(1.0f);

    return;
    }

  double dx;  //delta of position x: x_k+1 - x_k
  double dg;  //delta of gradient: g_k+1 - g_k
  double edg; //estimated delta of gradient: hessian_k * dx

  for (unsigned int i=0; i<numPara; i++)
    {
    dx = currentPositionRef[i] - this->m_PreviousPosition[i];
    dg = this->m_Gradient[i] - this->m_PreviousGradient[i];
    edg = m_LocalHessian[i] * dx;

    double dot_dg_dx = dg * dx;
    double dot_edg_dx = edg * dx;

    if (dot_dg_dx ==0 || dot_edg_dx == 0)
      {
      m_LocalHessian[i] = NumericTraits<double>::max();
      }
    else
      {
      double plus  = (dg * dg) / dot_dg_dx;
      double minus = (edg * edg) / dot_edg_dx;

      m_LocalHessian[i] = m_LocalHessian[i] + plus - minus;
      //m_LocalHessianInverse[i] = 1/newHessian;
      }
    }
}

/** Compute learning late from voxel shift*/
double QuasiNewtonObjectOptimizer
::EstimateLearningRate(ParametersType step)
{
  ParametersType parameters = this->GetCurrentPosition();

  ScalesType     scales = this->GetScales();

  double shift, learningRate;

  shift = this->m_Metric->ComputeMaximumVoxelShift(true, step);
  //shift = m_OptimizerHelper->ComputeMaximumVoxelShift(parameters, step);

  //initialize for the first time of executing EstimateLearningRate
  if (this->GetCurrentIteration() == 0 || m_MinimumVoxelShift == 0)
    {
    m_MinimumVoxelShift = shift * 1e-3;
    std::cout << " Initial learningRate = " << m_MaximumVoxelShift / shift << std::endl;
    }

  if (shift >= m_MinimumVoxelShift)
    {
    learningRate = m_MaximumVoxelShift / shift;
    }
  else //if the step yields almost zero voxel shift, it is not good to go
    {
    learningRate = 0;
    }
  return learningRate;
}

/** Debug the step sizes and turn angles */
void QuasiNewtonObjectOptimizer
::DebugStepSizeAndAngles(const char *debugLabel,
                         ParametersType lastStep,
                         ParametersType thisStep) const
{
  double product = inner_product(thisStep, lastStep);
  double lastNorm = vcl_sqrt(inner_product(lastStep, lastStep));
  double thisNorm = vcl_sqrt(inner_product(thisStep, thisStep));

  double cosine = product / lastNorm / thisNorm;
  double angle  = vcl_acos(cosine) / vnl_math::pi * 180;
  std::cout << debugLabel << " product = " << product
    << " angle = " << angle << " cosine = " << cosine
    << " lastNorm = " << lastNorm << " thisNorm = " << thisNorm
    << std::endl;

}

} // end namespace itk

#endif
