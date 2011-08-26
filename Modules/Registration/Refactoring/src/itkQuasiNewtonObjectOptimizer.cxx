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
//#ifndef __itkQuasiNewtonObjectOptimizer_hxx
//#define __itkQuasiNewtonObjectOptimizer_hxx

#include <iostream>
#include <iomanip>

#include "itkQuasiNewtonObjectOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkMacro.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"

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
  m_MinimumGradientNorm = 1e-20;
  m_MinimumValueChange = 1e-20;

  m_LocalHessian = NULL;
  m_LocalHessianInverse = NULL;

  m_LineSearchEnabled = true;
  m_OptimizerParameterEstimator = (OptimizerParameterEstimatorBase::Pointer)NULL;

  this->SetDebug(true);

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

  m_ValueAndDerivateEvaluated = false;

  if ( ! this->m_Metric->HasLocalSupport() )
    {
    if (m_OptimizerParameterEstimator.IsNotNull())
      {
      // initialize scales
      m_CurrentPosition = this->m_Metric->GetParameters();
      if (m_CurrentPosition[0] != m_CurrentPosition[0]) //checking NaN or #IND
        {
        itkExceptionMacro("QuasiNewtonObjectOptimizer: metric parameters are not defined.");
        }
      ScalesType scales(this->m_Metric->GetNumberOfParameters());

      m_OptimizerParameterEstimator->EstimateScales(scales);
      //m_CurrentPosition = this->m_Metric->GetParameters();
      this->SetScales(scales);
      std::cout << " Estimated scales = " << scales << std::endl;
      }
    }
  else
    {
    //no longer used for quasi-newton with localsupport
    //optimizer->SetLearningRate( 1.0 );
    //optimizer->SetScalarScale( 1.0 );
    //optimizer->SetUseScalarScale(true);
    }

  this->Superclass::StartOptimization();

}

/**
 * Resume optimization.
 */
void
QuasiNewtonObjectOptimizer
::ResumeOptimization()
{
  m_StopConditionDescription.str("");
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
  InvokeEvent( StartEvent() );

  m_Stop = false;
  while( ! m_Stop )
    {
    /* Compute value/derivative, using threader. */
    try
      {
      /* m_Gradient will be sized as needed by metric. If it's already
       * proper size, no new allocation is done. */
      if (!m_ValueAndDerivateEvaluated)
        {
        this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );

        }
      }
    catch ( ExceptionObject & err )
      {
      m_StopCondition = MetricError;
      m_StopConditionDescription << "Metric error";
      this->StopOptimization();

      // Pass exception to caller
      throw err;
      }

    /* Check if optimization has been stopped externally.
     * (Presumably this could happen from a multi-threaded client app?) */
    if ( m_Stop )
      {
      m_StopConditionDescription << "StopOptimization() called";
      break;
      }

    /* Advance one step along the gradient.
     * This will modify the gradient and update the transform. */
    this->AdvanceOneStep();

    /* Update and check iteration count */
    m_CurrentIteration++;
    if ( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                                 << ") exceeded.";
      m_StopCondition = MaximumNumberOfIterations;
      this->StopOptimization();
      break;
      }
    } //while (!m_Stop)
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
    if (m_OptimizerParameterEstimator.IsNull())
      {
      itkExceptionMacro("m_OptimizerParameterEstimator == NULL");
      }
    AdvanceOneLocalStep();
    return;
    }

  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  ScalesType scales = this->GetScales();

  double learningRate;
  bool   newtonStepException = false;

  this->m_CurrentPosition = this->m_Metric->GetParameters();

  if ( this->m_Gradient[0] != this->m_Gradient[0] ) //checking NaN
    {
    itkExceptionMacro("Gradient is undefined");
    }

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousPosition = this->m_CurrentPosition;
    m_PreviousGradient = this->m_Gradient;

    m_NewtonStep.SetSize(spaceDimension);
    m_Hessian.SetSize(spaceDimension, spaceDimension);
    m_HessianInverse.SetSize(spaceDimension, spaceDimension);
    }

  try
    {
    this->EstimateNewtonStep();
    }
  catch ( ExceptionObject & )
    {
    //This may happen with a singular hessian matrix
    std::cout << "Warning: exception in estimating Newton step." << std::endl;
    newtonStepException = true;
    }

  if (this->GetDebug())
    {
    std::cout << "m_CurrentPosition(:," << 1+this->GetCurrentIteration() << ") = " << this->m_CurrentPosition << "';" << std::endl;
    std::cout << "m_Value(" << 1+this->GetCurrentIteration() << ") = " << this->GetValue() << ";" << std::endl;
    std::cout << "m_NewtonStep(:," << 1+this->GetCurrentIteration() << ") = " << m_NewtonStep << "';" << std::endl;
    std::cout << "m_Gradient(:," << 1+this->GetCurrentIteration() << ") = " << m_Gradient << "';" << std::endl;
    }
  /** Save for the next iteration */
  m_PreviousPosition = this->GetCurrentPosition();
  m_PreviousGradient = this->GetGradient();

  /** If a Newton step is on the opposite direction of a gradient step, we'd
   * better use the gradient step. This happens when the second order
   * approximation produces a convex instead of an expected concave, or
   * vice versa.
   */
  if ( newtonStepException || inner_product(m_Gradient, m_NewtonStep) < 0 )
    {
    // using gradient step
    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small.";
      this->StopOptimization();
      return;
      }

    DerivativeType gradientStep = m_Gradient;
    for (unsigned int i=0; i<spaceDimension; i++)
      {
      gradientStep[i] = gradientStep[i] / scales[i];
      }
    learningRate = this->EstimateLearningRate(gradientStep);
    this->SetLearningRate( learningRate );

    if ( m_LineSearchEnabled )
      {
      //this->AdvanceWithSimpleLineSearch(gradientStep, learningRate);
      this->AdvanceWithStrongWolfeLineSearch(gradientStep, learningRate);
      this->m_ValueAndDerivateEvaluated = true;
      m_CurrentIteration--;
      }
    else
      {
      if (this->GetDebug())
        {
        std::cout << "using gradient, learningRate = " << learningRate << std::endl;
        }

      this->GradientDescentObjectOptimizer::AdvanceOneStep();
      this->m_ValueAndDerivateEvaluated = false;
      }
    }
  else
    {
    // Now a Newton step is on the consistent direction of a gradient step
    // using Newton step
    double newtonNorm = m_NewtonStep.two_norm();
    if (newtonNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the Newton step is too small.";
      this->StopOptimization();
      return;
      }

    learningRate = this->EstimateLearningRate(m_NewtonStep);
    learningRate = vnl_math_min(learningRate, 1.0);
    this->SetLearningRate( learningRate );

    if ( m_LineSearchEnabled )
      {
      //this->AdvanceWithSimpleLineSearch(m_NewtonStep, learningRate);
      this->AdvanceWithStrongWolfeLineSearch(m_NewtonStep, learningRate);
      this->m_ValueAndDerivateEvaluated = true;
      m_CurrentIteration--;
      }
    else
      {
      if (this->GetDebug())
        {
        std::cout << "using newton, learningRate = " << learningRate << std::endl;
        }

      DerivativeType step(spaceDimension);
      for ( unsigned int j = 0; j < spaceDimension; j++ )
        {
        step[j] = this->m_LearningRate * this->m_NewtonStep[j];
        }
      this->m_Metric->UpdateTransformParameters( step );
      this->m_ValueAndDerivateEvaluated = false;

      this->InvokeEvent( IterationEvent() );
      } // without line search
    } // using newton step
}

/*************************************************
 * Do backtracking line search on the Newton direction
 *************************************************/
void QuasiNewtonObjectOptimizer
::AdvanceWithSimpleLineSearch(ParametersType direction, double maxStepSize)
{
  double optimalDirection = -1.0; //maximizing

  double t = 1.0, beta = 0.75;
  double c1 = 1e-4;

  double oldValue = this->GetValue();
  double stepChange = inner_product(direction, m_Gradient);

  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  ParametersType step(spaceDimension);

  ParametersType initPosition = this->m_Metric->GetParameters();
  ParametersType tempPosition(spaceDimension);

  while (true)
    {
    tempPosition = this->m_Metric->GetParameters();
    step = initPosition + direction * t - tempPosition;

    this->m_Metric->UpdateTransformParameters( step );
    this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
    m_CurrentIteration++;

    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small in line search.";
      this->StopOptimization();
      return;
      }

    if ((this->m_Value - (oldValue + c1 * t * stepChange)) * optimalDirection <= 0)
      {
      break;
      }
    /* Update and check iteration count */
    if ( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                                 << ") exceeded.";
      m_StopCondition = MaximumNumberOfIterations;
      this->StopOptimization();
      return;
      }

    t *= beta;

    }

  //std::cout << "line search step size = " << t << std::endl;

  this->InvokeEvent( IterationEvent() );

}

/*************************************************
 * Do backtracking line search on the Newton direction
 *************************************************/
double QuasiNewtonObjectOptimizer
::AdvanceWithStrongWolfeLineSearch(ParametersType direction, double maxStepSize)
{
  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  double optimalDirection = -1.0; //maximizing

  double tmax = maxStepSize, t0 = 0, topt;
  double t1 = t0, t2 = tmax / 2.0;

  double c1 = 1e-4, c2 = 0.9;

  double f0, f1, f2;
  double g0, g1, g2; //derivative w.r.t the step size t

  f0 = optimalDirection * this->GetValue();
  g0 = optimalDirection * inner_product(direction, this->m_Gradient);

  f1 = f0;
  g1 = g0;

  int loop = 1;
  ParametersType initPosition = this->m_Metric->GetParameters();
  ParametersType tempPosition(spaceDimension);
  ParametersType deltaPosition(spaceDimension);

  while (true)
    {
    tempPosition = this->m_Metric->GetParameters();
    deltaPosition = initPosition + t2 * direction - tempPosition;

    this->m_Metric->UpdateTransformParameters( deltaPosition );
    this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
    m_CurrentIteration++;

    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small in line search.";
      this->StopOptimization();

      topt = t2;
      return topt;
      }

    f2 = optimalDirection * this->m_Value;
    g2 = optimalDirection * inner_product(direction, this->m_Gradient);

    if (f2 > f0 + c1 * t2 * g0 || ( loop > 1 && f2 >= f1 ))
      {
      topt = this->LineSearchZoom(initPosition, f0, g0, direction, t1, t2);
      return topt;
      }

    if (vcl_abs(g2) <= -c2 * g0)
      {
      topt = t2;
      return topt;
      }

    if (g2 >= 0)
      {
      topt = this->LineSearchZoom(initPosition, f0, g0, direction, t2, t1);
      this->InvokeEvent( IterationEvent() );
      return topt;
      }
    t1 = t2;
    f1 = f2;
    g1 = g2;
    t2 = t2 + (tmax - t2) / 2;

    loop++;

    /* Update and check iteration count */
    if ( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                                 << ") exceeded in line search. ";
      m_StopCondition = MaximumNumberOfIterations;
      this->StopOptimization();
      topt = t2;
      return topt;
      }
    } //while

}

double QuasiNewtonObjectOptimizer
::LineSearchZoom(ParametersType initPosition, double f0, double g0, ParametersType direction, double tlow, double thigh)
{
  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  double optimalDirection = -1.0; //maximizing

  double c1 = 1e-4, c2 = 0.9;
  double t2, topt;
  double f2, g2;
  double flow, glow;

  double tempValue, oldValue = f0;

  ParametersType tempPosition(spaceDimension);
  ParametersType deltaPosition(spaceDimension);
  DerivativeType tempGradient(spaceDimension);

  int loop = 0;
  while (true)
    {
    loop++;
    t2 = (tlow + thigh) / 2.0;

    tempPosition = this->m_Metric->GetParameters();
    deltaPosition = initPosition + t2 * direction - tempPosition;

    this->m_Metric->UpdateTransformParameters( deltaPosition );
    this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
    m_CurrentIteration++;

    if (this->GetDebug())
      {
      std::cout << "LineSearchZoom: position=" << initPosition + t2 * direction << std::endl;
      std::cout << "LineSearchZoom: m_Value=" << this->m_Value << std::endl;
      std::cout << "LineSearchZoom: m_Gradient=" << this->m_Gradient << std::endl;
      }

    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small in Wolfe line search.";
      this->StopOptimization();

      topt = t2;
      return topt;
      }

    f2 = optimalDirection * this->m_Value;
    g2 = optimalDirection * inner_product(direction, this->m_Gradient);

    if (loop >= 2 && vcl_abs(oldValue - f2) < m_MinimumValueChange)
      {
      m_StopCondition = StepTooSmall;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the value change is too small in Wolfe line search.";
      this->StopOptimization();

      topt = t2;
      return topt;
      }

    oldValue = f2;

    if ( f2 > f0 + c1 * t2 * g0 )
      {
      thigh = t2;
      }
    else
      {
      tempPosition = this->m_Metric->GetParameters();
      deltaPosition = initPosition + tlow * direction - tempPosition;

      this->m_Metric->UpdateTransformParameters( deltaPosition );
      this->m_Metric->GetValueAndDerivative( tempValue, tempGradient );
      m_CurrentIteration++;

      flow = optimalDirection * tempValue;
      glow = optimalDirection * inner_product(direction, tempGradient);

      if ( f2 >= flow ) // f2 == flow when |thigh - tlow| is small?
        {
        thigh = t2;
        }
      else
        {
        if ( vcl_abs(g2) <= -c2 * g0 )
          {
          //reset the parameters to initPosition + t2 * direction
          //this->m_Value and this->m_Gradient are already evaluated
          tempPosition = this->m_Metric->GetParameters();
          deltaPosition = initPosition + t2 * direction - tempPosition;
          this->m_Metric->UpdateTransformParameters( deltaPosition );

          this->InvokeEvent( IterationEvent() );

          topt = t2;
          return topt;
          }

        if (g2 * (thigh - tlow) >= 0)
          {
          thigh = tlow;
          }
        tlow = t2;
        }
      }

    if ( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                               << ") exceeded in Wolfe line search. ";
      m_StopCondition = MaximumNumberOfIterations;
      this->StopOptimization();

      topt = t2;

      //reset the parameters to initPosition + t2 * direction
      //this->m_Value and this->m_Gradient are already evaluated
      tempPosition = this->m_Metric->GetParameters();
      deltaPosition = initPosition + t2 * direction - tempPosition;
      this->m_Metric->UpdateTransformParameters( deltaPosition );

      return topt;
      }
    } //while
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
  const unsigned int imageDimension = 2;
  const unsigned int imageSize = spaceDimension / imageDimension;
  ScalesType scales = this->GetScales();

  bool   newtonStepException = false;

  // Use reference to save memory copy
  const ParametersType & currentPositionRef = this->m_Metric->GetParameters();

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousPosition = currentPositionRef;
    m_PreviousGradient = this->m_Gradient;

    m_NewtonStep.SetSize(spaceDimension);

    m_LocalHessian = new LocalHessianType[imageSize];
    m_LocalHessianInverse = new LocalHessianType[imageSize];

    unsigned int imgDim = this->m_OptimizerParameterEstimator->GetImageDimension();
    for (unsigned int i=0; i<imageSize; i++)
      {
      m_LocalHessian[i].SetSize(imgDim, imgDim);
      m_LocalHessianInverse[i].SetSize(imgDim, imgDim);
      }
    }

  try
    {
    this->EstimateLocalNewtonStep();
    }
  catch ( ExceptionObject )
    {
    //This may happen with a singular hessian matrix
    std::cout << "Warning: exception in estimating Newton step." << std::endl;
    newtonStepException = true;
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

  double gradientLearningRate = 100*m_MaximumVoxelShift / maxGradient;
  double newtonLearningRate   = 100*m_MaximumVoxelShift / maxNewtonStep;
  newtonLearningRate = vnl_math_min(newtonLearningRate, 1.0);

  if (this->GetDebug())
    {
    std::cout << "Iteration = " << this->GetCurrentIteration() << std::endl;
    std::cout << "spaceDimension = " << spaceDimension << std::endl;
    std::cout << "newtonLearningRate = " << newtonLearningRate << std::endl;
    std::cout << "gradientLearningRate = " << gradientLearningRate << std::endl;
    }

  for ( unsigned int i=0; i<imageSize; i++ )
    {
    /** If a Newton step is on the opposite direction of a gradient step, we'd
     * better use the gradient step. This happens when the second order
     * approximation produces a convex instead of an expected concave, or
     * vice versa.
     */
    double dotProduct = 0;
    for (unsigned int d=0; d<imageDimension; d++)
      {
      dotProduct += m_Gradient[i*imageDimension + d] * m_NewtonStep[i*imageDimension + d];
      }
    if ( dotProduct <= 0 )
      {
      for (unsigned int d=0; d<imageDimension; d++)
        {
        m_NewtonStep[i*imageDimension + d] = m_Gradient[i*imageDimension + d] * gradientLearningRate;
        }
      }
    else
      {
      for (unsigned int d=0; d<imageDimension; d++)
        {
        m_NewtonStep[i*imageDimension + d] = m_NewtonStep[i*imageDimension + d] * newtonLearningRate;
        }
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

  if (this->GetCurrentIteration() == 0 ||
      m_Hessian[0][0] == NumericTraits<double>::max() ||
      m_HessianInverse[0][0] == NumericTraits<double>::max())
    {
    m_NewtonStep = m_Gradient; //use gradient step

    // Initialize Hessian to identity matrix
    m_Hessian.Fill(0.0f);
    m_HessianInverse.Fill(0.0f);
    for (unsigned int i=0; i<numPara; i++)
      {
      m_Hessian[i][i] = 1.0; //identity matrix
      m_HessianInverse[i][i] = 1.0; //identity matrix
      }
    }
  else
    {
    sdx = m_HessianInverse * m_Gradient;
    //this->m_NewtonStep = -1.0 * sdx; //minimize
    this->m_NewtonStep = sdx; //maximize
    }
}

/** Estimate Hessian matrix */
void QuasiNewtonObjectOptimizer
::EstimateHessian()
{
  if (this->GetCurrentIteration() == 0)
    {
    return;
    }

  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  ParametersType dx(numPara);  //delta of position x: x_k+1 - x_k
  ParametersType dg(numPara);  //delta of gradient: g_k+1 - g_k
  ParametersType edg(numPara); //estimated delta of gradient: hessian_k * dx

  dx = this->m_CurrentPosition - this->m_PreviousPosition;
  //dg = this->m_Gradient - this->m_PreviousGradient; //minimize
  dg = this->m_PreviousGradient - this->m_Gradient; //maximize
  edg = m_Hessian * dx;

  double dot_dg_dx = inner_product(dg, dx);
  double dot_edg_dx = inner_product(edg, dx);

  if (dot_dg_dx ==0 || dot_edg_dx == 0)
    {
    m_Hessian[0][0] = NumericTraits<double>::max();
    return;
    }

  vnl_matrix<double> plus  = outer_product(dg, dg) / dot_dg_dx;
  vnl_matrix<double> minus = outer_product(edg, edg) / dot_edg_dx;
  vnl_matrix<double> newHessian = m_Hessian + plus - minus;

  m_Hessian         = newHessian;

  if ( vcl_abs(vnl_determinant(newHessian)) == 0 )
    {
    m_HessianInverse[0][0] = NumericTraits<double>::max();
    }
  else
    {
    m_HessianInverse = vnl_matrix_inverse<double>(newHessian);
    }
}

/** Estimate Hessian step */
void QuasiNewtonObjectOptimizer
::EstimateLocalNewtonStep()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();
  const unsigned int imageDimension = 2;
  const unsigned int imageSize = numPara / imageDimension;

  ParametersType localGradient(imageDimension);
  ParametersType localNewtonStep(imageDimension);

  // Estimate Hessian
  EstimateLocalHessian();

  // Compute the Newton step
  for (unsigned int i=0; i<imageSize; i++)
    {
    if (this->GetCurrentIteration() == 0 ||
        m_LocalHessian[i][0][0] == NumericTraits<double>::max() ||
        m_LocalHessianInverse[i][0][0] == NumericTraits<double>::max())
      {
      m_LocalHessian[i].Fill(0.0);
      m_LocalHessianInverse[i].Fill(0.0);
      for (unsigned int d=0; d<imageDimension; d++)
        {
        m_NewtonStep[i*imageDimension + d] = m_Gradient[i*imageDimension + d]; //use gradient step
        m_LocalHessian[i][d][d] = 1; //reset to identity
        m_LocalHessianInverse[i][d][d] = 1; //reset to identity
        }
      }
    else
      {
      for (unsigned int d=0; d<imageDimension; d++)
        {
        localGradient[d] = m_Gradient[i*imageDimension + d];
        }
      //m_NewtonStep[i] = - m_Gradient[i] / m_LocalHessian[i];
      localNewtonStep = m_LocalHessianInverse[i] * localGradient;

      for (unsigned int d=0; d<imageDimension; d++)
        {
        m_NewtonStep[i*imageDimension + d] = localNewtonStep[d];
        }
      }
    }

}

/** Estimate Hessian matrix */
void QuasiNewtonObjectOptimizer
::EstimateLocalHessian()
{
  if (this->GetCurrentIteration() == 0)
    {
    return;
    }

  unsigned int numPara = this->m_Metric->GetNumberOfParameters();
  const unsigned int imageDimension = 2;
  const unsigned int imageSize = numPara / imageDimension;

  // Use reference to save memory copy
  const ParametersType & currentPositionRef = this->m_Metric->GetParameters();

  ParametersType dx(imageDimension);  //delta of position x: x_k+1 - x_k
  ParametersType dg(imageDimension);  //delta of gradient: g_k+1 - g_k
  ParametersType edg(imageDimension); //estimated delta of gradient: hessian_k * dx

  for (unsigned int i=0; i<imageSize; i++)
    {
    for (unsigned int j=0; j<imageDimension; j++)
      {
      dx[j] = currentPositionRef[i*imageDimension+j] - this->m_PreviousPosition[i*imageDimension+j];
      dg[j] = this->m_PreviousGradient[i*imageDimension+j] - this->m_Gradient[i*imageDimension+j]; //maximize
      }
    edg = m_LocalHessian[i] * dx;

    double dot_dg_dx = inner_product(dg, dx);
    double dot_edg_dx = inner_product(edg, dx);

    if (dot_dg_dx ==0 || dot_edg_dx == 0)
      {
      m_LocalHessian[i][0][0] = NumericTraits<double>::max();
      }
    else
      {
      vnl_matrix<double> plus  = outer_product(dg, dg) / dot_dg_dx;
      vnl_matrix<double> minus = outer_product(edg, edg) / dot_edg_dx;
      vnl_matrix<double> newHessian = m_LocalHessian[i] + plus - minus;

      m_LocalHessian[i] = newHessian;

      if ( vcl_abs(vnl_determinant(newHessian)) <= 0.5 )
        {
        m_LocalHessianInverse[i][0][0] = NumericTraits<double>::max();
        }
      else
        {
        m_LocalHessianInverse[i] = vnl_matrix_inverse<double>(newHessian);
        }
      }
    }
}

/** Compute learning late from voxel shift*/
double QuasiNewtonObjectOptimizer
::EstimateLearningRate(ParametersType step)
{
  if (m_OptimizerParameterEstimator.IsNull())
    {
    return 1;
    }

  ParametersType parameters = this->GetCurrentPosition();

  ScalesType     scales = this->GetScales();

  double shift, learningRate;

  shift = m_OptimizerParameterEstimator->ComputeMaximumVoxelShift(step);

  if (this->GetCurrentIteration() == 0)
    {
    std::cout << " Initial learningRate = " << m_MaximumVoxelShift / shift << std::endl;
    }

  if (shift > 0)
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

//#endif
