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
#ifndef _itkQuasiNewtonOptimizer_txx
#define _itkQuasiNewtonOptimizer_txx

#include "itkQuasiNewtonOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkMacro.h"

namespace itk
{
/**
 * Constructor
 */
QuasiNewtonOptimizer
::QuasiNewtonOptimizer()
{
  itkDebugMacro("Constructor");

  m_LineSearchEnabled = false;

}

void
QuasiNewtonOptimizer
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * Advance one Step following the gradient direction
 */
void
QuasiNewtonOptimizer
::AdvanceOneStep(void)
{
  itkDebugMacro("AdvanceOneStep");

  double direction;
  if ( this->m_Maximize )
    {
    direction = 1.0;
    }
  else
    {
    direction = -1.0;
    }

  const unsigned int spaceDimension =  m_CostFunction->GetNumberOfParameters();
  unsigned int curIt = this->GetCurrentIteration();
  double learningRate;

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousPosition = m_CurrentPosition;
    m_PreviousGradient = m_Gradient;
    }

  try
    {
    this->EstimateNewtonStep();
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

  if (this->GetDebug() && this->GetCurrentIteration() > 0)
    {
    this->DebugStepSizeAndAngles("PreviousGradient Gradient",
                                  m_PreviousGradient, m_Gradient);
    this->DebugStepSizeAndAngles("PreviousNewtonStep NewtonStep",
                                  m_PreviousNewtonStep, m_NewtonStep);
    this->DebugStepSizeAndAngles("Gradient NewtonStep",
                                  m_Gradient, m_NewtonStep);
    }

  /** Save for the next iteration */
  m_PreviousPosition = this->GetCurrentPosition();
  m_PreviousGradient = this->GetGradient();
  m_PreviousNewtonStep = m_NewtonStep;

  /** If a Newton step is on the opposite direction of a gradient step, we'd
   * better use the gradient step. This happens when the second order
   * approximation produces a convex instead of an expected concave, or
   * vice versa.
   */
  if ( (!this->m_Maximize && inner_product(m_Gradient, m_NewtonStep) >= 0) ||
       ( this->m_Maximize && inner_product(m_Gradient, m_NewtonStep) <= 0) )
    {
    ParametersType step = m_Gradient;
    step = direction * step;
    learningRate = this->Superclass::EstimateLearningRate(step);
    if ( this->GetDebug() )
      {
      std::cout << curIt << " UseGradient " << std::endl;
      std::cout << " Using gradient step" << std::endl;
      std::cout << " learningRate" << curIt << " = " << learningRate << std::endl;
      std::cout << " m_CurrentPosition" << curIt << " = " << m_CurrentPosition
        << std::endl;
      std::cout << " m_Gradient"   << curIt << " = " << m_Gradient << std::endl;
      }
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
    this->GradientDescentOptimizer::AdvanceOneStep();
    return;
    }

  // Now a Newton step is on the consistent direction of a gradient step
  learningRate = this->Superclass::EstimateLearningRate(m_NewtonStep);
  learningRate = vnl_math_min(learningRate, 1.0);

  if ( this->GetDebug() )
    {
    std::cout << curIt << " UseNewton   " << std::endl;
    std::cout << " m_CurrentPosition" << curIt << " = " << m_CurrentPosition
      << std::endl;
    std::cout << " m_Gradient"   << curIt << " = " << m_Gradient << std::endl;
    std::cout << " m_NewtonStep" << curIt << " = " << m_NewtonStep << std::endl;
    std::cout << " m_ScaledNewtonStep" << curIt << " = " << m_ScaledNewtonStep
      << std::endl;
    std::cout << " learningRate" << curIt << " = " << learningRate << std::endl;

    ParametersType scaledPos(spaceDimension);
    DerivativeType scaledGra(spaceDimension);
    this->ScalePosition(m_CurrentPosition, scaledPos);
    this->ScaleDerivative(m_Gradient, scaledGra);
    std::cout << " scaledPos" << curIt << " = " << scaledPos << std::endl;
    std::cout << " scaledGra" << curIt << " = " << scaledGra << std::endl;
    }

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

  if ( ! m_LineSearchEnabled )
    {
    DerivativeType curGradient = m_Gradient;
    m_Gradient = -m_NewtonStep;
    this->GradientDescentOptimizer::AdvanceOneStep();
    m_Gradient = curGradient;
    }
  else
    {
    this->LineSearch();
    }

}
  /*************************************************
   * Do backtracking line search on the Newton direction
   *************************************************/
void QuasiNewtonOptimizer::LineSearch()
{
  double direction = m_Maximize ? 1.0 : -1.0;

  double learningRate = m_LearningRate;
  const unsigned int spaceDimension =  m_CostFunction->GetNumberOfParameters();
  unsigned int curIt = this->GetCurrentIteration();

  const ParametersType & currentPosition = this->GetCurrentPosition();
  ScalesType scales = this->GetScales();

  ParametersType scaledPosition(spaceDimension);
  ParametersType scaledCurrentPosition(spaceDimension);
  ParametersType newPosition(spaceDimension);
  ParametersType scaledStep(spaceDimension);

  DerivativeType scaledGradient(spaceDimension);
  DerivativeType newGradient(spaceDimension);
  DerivativeType scaledNewGradient(spaceDimension);

  for ( unsigned int j = 0; j < spaceDimension; j++ )
    {
    scaledStep[j] = learningRate * m_ScaledNewtonStep[j];
    }

  this->ScalePosition(m_CurrentPosition, scaledCurrentPosition);
  this->ScaleDerivative(m_Gradient, scaledGradient);

  double newValue, oldValue;
  oldValue = this->GetValue();

  double t = 1.0, beta = 0.75;
  double c1 = 1e-4, c2 = 0.9;
  double stepChange = inner_product(scaledStep, scaledGradient);

  for (int searchCount = 0; searchCount < 20; searchCount++)
    {
    for ( unsigned int j = 0; j < spaceDimension; j++ )
      {
      scaledPosition[j] = scaledCurrentPosition[j] + t * scaledStep[j];
      }
    this->ScaleBackPosition(scaledPosition, newPosition);
    m_CostFunction->GetValueAndDerivative(newPosition, newValue, newGradient);

    //Wolfe condition I
    if ((newValue - (oldValue + c1 * t * stepChange)) * direction >= 0)
      {
      break;
      }
    /*
    //Wolfe condition II
    this->ScaleDerivative(newGradient, scaledNewGradient);
    if ( (inner_product(scaledStep, scaledNewGradient) - c2 * stepChange)
         * direction <= 0)
      {
      break;
      }
    */
    t *= beta;
    }

  this->SetCurrentPosition(newPosition);

  this->InvokeEvent( IterationEvent() );

}
/** Estimate Hessian step */
void QuasiNewtonOptimizer
::EstimateNewtonStep()
{
  int numPara = m_CurrentPosition.size();
  ParametersType sx1(numPara);
  ParametersType sx2(numPara);
  DerivativeType sg1(numPara);
  DerivativeType sg2(numPara);

  // Translate to the scaled space
  this->ScalePosition(m_PreviousPosition, sx1);
  this->ScalePosition(m_CurrentPosition,  sx2);
  this->ScaleDerivative(m_PreviousGradient, sg1);
  this->ScaleDerivative(m_Gradient,         sg2);

  // Estimate Hessian in the scaled space
  EstimateHessian(sx1, sx2, sg1, sg2);

  // Compute the Newton step
  ParametersType sdx(numPara);
  sdx = m_HessianInverse * sg2;
  sdx = -1.0 * sdx;
  m_ScaledNewtonStep = sdx;

  // Translate the step back into the original space
  ParametersType dx(numPara);
  this->ScaleBackDerivative(sdx, dx);
  m_NewtonStep = dx;

}

/** Estimate Hessian matrix */
void QuasiNewtonOptimizer
::EstimateHessian(ParametersType x1, ParametersType x2,
                  DerivativeType g1, DerivativeType g2)
{
  int numPara = x1.size();

  // Initialize Hessian to identity matrix
  if ( this->GetCurrentIteration() == 0 )
    {
    m_Hessian.SetSize(numPara, numPara);
    m_Hessian.Fill(0.0f);
    m_HessianInverse.SetSize(numPara, numPara);
    m_HessianInverse.Fill(0.0f);

    for (int i=0; i<numPara; i++)
      {
      m_Hessian[i][i] = 1.0; //identity matrix
      m_HessianInverse[i][i] = 1.0; //identity matrix
      }
    return;
    }

  HessianType oldHessian(numPara, numPara);
  oldHessian = m_Hessian;

  ParametersType dx(numPara);  //delta of position x: x_k+1 - x_k
  ParametersType dg(numPara);  //delta of gradient: g_k+1 - g_k
  ParametersType edg(numPara); //estimated delta of gradient: hessian_k * dx

  dx = x2 - x1;
  dg = g2 - g1;
  edg = oldHessian * dx;

  double dot_dg_dx = inner_product(dg, dx);
  double dot_edg_dx = inner_product(edg, dx);

  if (dot_dg_dx ==0 || dot_edg_dx == 0)
    {
    itkExceptionMacro(<< "Division by zero in Quasi-Newton step. ");
    return;
    }

  vnl_matrix<double> plus  = outer_product(dg, dg) / dot_dg_dx;
  vnl_matrix<double> minus = outer_product(edg, edg) / dot_edg_dx;
  vnl_matrix<double> newHessian = oldHessian + plus - minus;


  if ( this->GetDebug() )
    {
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "g1 = " << g1 << std::endl;
    std::cout << "g2 = " << g2 << std::endl;
    std::cout << "oldHessian = [" << (vnl_matrix<double>)oldHessian << "]"
      << std::endl;
    std::cout << "newHessian = [" << newHessian << "]" << std::endl;
    }

  m_Hessian         = newHessian;
  m_HessianInverse  = vnl_matrix_inverse<double>(newHessian);

}

/** Translate the parameters into the scaled space */
void QuasiNewtonOptimizer::ScalePosition(ParametersType p1, ParametersType &p2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<p1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      p2[i] = p1[i] * scale;
    }
}

/** Translate the parameters into the original space */
void QuasiNewtonOptimizer::ScaleBackPosition(ParametersType p1, ParametersType &p2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<p1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      p2[i] = p1[i] / scale;
    }
}

/** Translate the derivative into the scaled space */
void QuasiNewtonOptimizer::ScaleDerivative(DerivativeType g1, DerivativeType &g2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<g1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      g2[i] = g1[i] / scale;
    }
}


/** Translate the derivative into the original space */
void QuasiNewtonOptimizer::ScaleBackDerivative(DerivativeType g1, DerivativeType &g2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<g1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      g2[i] = g1[i] * scale;
    }
}

} // end namespace itk

#endif
