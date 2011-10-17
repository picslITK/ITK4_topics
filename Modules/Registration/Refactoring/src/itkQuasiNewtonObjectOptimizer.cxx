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
#include <vnl/algo/vnl_determinant.h>

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
  m_MaximumNewtonVoxelShift = 5.0;

  m_MinimumGradientNorm = 1e-15;
  m_MinimumValueChange = 1e-15;

  m_LineSearchEnabled = false;
  m_AlgorithmOption = "qn"; //quasi-newton, others may be lqn cg gd
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
  if (m_AlgorithmOption.empty())
    {
    m_AlgorithmOption = "qn";
    }
  m_ValueAndDerivateEvaluated = false;

  if ( ! this->m_Metric->HasLocalSupport() )
    {
    if (m_OptimizerParameterEstimator.IsNotNull())
      {
      // initialize scales
      ScalesType scales(this->m_Metric->GetNumberOfLocalParameters());

      m_OptimizerParameterEstimator->EstimateScales(scales);
      //for (int s=0; s<4; s++) scales[s] = 9801;
      //for (int s=4; s<6; s++) scales[s] = 1;

      this->SetScales(scales);
      std::cout << " Estimated scales = " << scales << std::endl;
      }
    }
  else
    {
    itkExceptionMacro("To have local support, please use QuasiNewtonLocalSupportObjectOptimizer instead.");
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

  if (m_AlgorithmOption.compare("cg") == 0)
    {
    this->AdvanceWithConjugateGradient();
    return;
    }
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
        m_PastPositions.push_back(this->m_Metric->GetParameters());
        m_PastValues.push_back(this->m_Value);
        m_PastGradients.push_back(this->m_Gradient);
        }
      }
    catch ( ExceptionObject & err )
      {
      m_StopCondition = COSTFUNCTION_ERROR;
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
      m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
      //this->AdvanceWithNeighborhood();
      //this->AdvanceWithConjugateGradient();
      this->StopOptimization();
      break;
      }
    } //while (!m_Stop)
}

void
QuasiNewtonObjectOptimizer
::AdvanceOneStep(void)
{
  itkDebugMacro("AdvanceOneStep");

  if ( this->m_Metric->HasLocalSupport() )
    {
    itkExceptionMacro("To have local support, please use QuasiNewtonLocalSupportObjectOptimizer instead.");
    return;
    }

  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  ScalesType scales = this->GetScales();

  double learningRate;
  bool   newtonStepException = false;

  this->m_CurrentPosition = this->m_Metric->GetParameters();

  if (this->GetDebug())
    {
    m_HistoryValues.push_back(m_Value);
    m_HistoryPositions.push_back(m_CurrentPosition);
    }

  if ( this->m_Gradient[0] != this->m_Gradient[0] ) //checking NaN
    {
    itkExceptionMacro("Gradient is undefined");
    }

  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousValue = this->m_Value;
    m_PreviousPosition = this->m_CurrentPosition;
    m_PreviousGradient = this->m_Gradient;

    m_BestValue = this->m_Value;
    m_BestPosition = this->m_CurrentPosition;
    m_BestIteration = this->GetCurrentIteration();

    m_NewtonStep.SetSize(spaceDimension);
    m_Hessian.SetSize(spaceDimension, spaceDimension);
    m_HessianInverse.SetSize(spaceDimension, spaceDimension);
    }

  if (m_BestValue < this->m_Value)
    {
    m_BestValue = this->m_Value;
    m_BestPosition = this->m_CurrentPosition;
    m_BestIteration = this->GetCurrentIteration();
    }
  else
    {
    if ( this->GetCurrentIteration() >= m_BestIteration + 50)
      {
      ParametersType backStep;
      backStep = m_BestPosition - this->m_Metric->GetParameters();
      this->m_Metric->UpdateTransformParameters( backStep );
      this->m_CurrentPosition = m_BestPosition;
      this->m_Value = m_BestValue;

      m_StopCondition = STEP_TOO_SMALL;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the value change is too small.";
      this->StopOptimization();
      return;
      }
    }

  if ( this->GetCurrentIteration() > 0
    && vcl_abs(m_PreviousValue - m_Value) < m_MinimumValueChange )
    {
    m_StopCondition = STEP_TOO_SMALL;
    m_StopConditionDescription << "Optimization stops after "
                               << this->GetCurrentIteration()
                               << " iterations since"
                               << " the value change is too small.";
    this->StopOptimization();
    return;
    }

  try
    {
    //Estimate Quasi-Newton step
    if (m_AlgorithmOption.compare("lqn") == 0)
      {
      this->EstimateNewtonStepWithLBFGS();
      }
    else
      {
      this->EstimateNewtonStep();
      }
    }
  catch ( ExceptionObject & )
    {
    //This may happen with a singular hessian matrix
    std::cout << "Warning: exception in estimating Newton step." << std::endl;
    newtonStepException = true;
    }

  if (this->GetDebug())
    {
    int iter = 1 + this->GetCurrentIteration();
    std::cout << "m_CurrentPosition(" << iter << ",:) = " << this->m_CurrentPosition << "';" << std::endl;
    std::cout << "m_Value(" << iter << ") = " << this->GetValue() << ";" << std::endl;
    std::cout << "m_NewtonStep(" << iter << ",:) = " << m_NewtonStep << "';" << std::endl;
    std::cout << "m_Gradient(" << iter << ",:) = " << m_Gradient << "';" << std::endl;
    }
  /** Save for the next iteration */
  m_PreviousValue = this->GetValue();
  m_PreviousPosition = this->m_CurrentPosition;
  m_PreviousGradient = this->GetGradient();

  /** If a Newton step is on the opposite direction of a gradient step, we'd
   * better use the gradient step. This happens when the second order
   * approximation produces a convex instead of an expected concave, or
   * vice versa.
   */
  double gradientNewtonProduct = inner_product(m_Gradient, m_NewtonStep);
  double gradientNorm = m_Gradient.two_norm();
  double newtonNorm = m_NewtonStep.two_norm();
  //double cosine = gradientNewtonProduct / gradientNorm / newtonNorm;
  bool gradientDescentOnly = (m_AlgorithmOption.compare("gd") == 0);

  if ( gradientDescentOnly || newtonStepException || gradientNewtonProduct < 0 )
    {
    // using gradient step
    this->ResetNewtonStep();

    //double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = STEP_TOO_SMALL;
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
      this->AdvanceWithStrongWolfeLineSearch(gradientStep, learningRate*10, learningRate);
      this->m_ValueAndDerivateEvaluated = true;

      //m_CurrentIteration was added in line search for parameter changes.
      //it is decreased here and will be added in the end of ResumeOptimization.
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
    //double newtonNorm = m_NewtonStep.two_norm();
    if (newtonNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = STEP_TOO_SMALL;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the Newton step is too small.";
      this->StopOptimization();
      return;
      }

    learningRate = this->EstimateLearningRate(m_NewtonStep, m_MaximumNewtonVoxelShift);
    learningRate = vnl_math_min(learningRate, 1.0);
    this->SetLearningRate( learningRate );

    if ( m_LineSearchEnabled )
      {
      learningRate = 1.0;
      this->AdvanceWithStrongWolfeLineSearch(m_NewtonStep, learningRate*10, learningRate);
      this->m_ValueAndDerivateEvaluated = true;

      //m_CurrentIteration was added in line search for parameter changes.
      //it is decreased here and will be added in the end of ResumeOptimization.
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

/**
 * Advance one Step: if the Quasi-Newton direction is consistent
 * with the gradient direction, follow the Quasi-Newton direction,
 * otherwise, follow the gradient direction.
 */
void
QuasiNewtonObjectOptimizer
::AdvanceWithNeighborhood(void)
{
  itkDebugMacro("AdvanceDebugStep");
  const unsigned int pnum =  this->m_Metric->GetNumberOfParameters();
  const unsigned int rows = m_HistoryPositions.size();
  itk::Array2D<double> ypos(rows, pnum);
  itk::Array2D<double> xpos(rows, pnum);
  itk::Array<double>   ypostmp(pnum);
  itk::Array<double>   xpostmp(pnum);
  itk::Array<double>   xmean(pnum);
  itk::Array2D<double> coeff(pnum, pnum);
  itk::Array2D<double> coeffBack(pnum, pnum);

  xmean.fill(0.0);
  for (UnInt i=0; i<rows; i++)
    {
    xmean += m_HistoryPositions[i];
    }
  xmean /= rows;

  for (UnInt i=0; i<rows; i++)
    {
    xpos.set_row(i, m_HistoryPositions[i]-xmean);
    }

  vnl_svd<double> xsvd(xpos);
  std::cout << "DebugXMean = " << xmean << ";" << std::endl;
  for (UnInt i=0; i<pnum; i++)
    {
    //make "grep DebugSvdV" work
    std::cout << "DebugSvdV(" << i+1 << ",:) = [" << xsvd.V().get_row(i) << "];" << std::endl;
    }
  std::cout << "DebugSvdW = " << xsvd.W() << ";" << std::endl;
  coeff = xsvd.V();
  coeffBack = coeff.transpose();
  for (UnInt i=0; i<rows; i++)
    {
    ypos.set_row(i, (xpos.get_row(i)) * coeff);
    }

  int dimx = 201, dimy = 201;
  int cx = 0, cy = 1;
  double smallUnit = 2; //4;
  double yposx_low  = ypos[rows-1][cx] - smallUnit*2;
  double yposx_high = ypos[rows-1][cx] + smallUnit*2;
  double yposy_low  = ypos[rows-1][cy] - smallUnit;
  double yposy_high = ypos[rows-1][cy] + smallUnit;
  yposx_low  = ypos[1-1][cx] - smallUnit*2;
  yposx_high = ypos[rows-1][cx] + smallUnit*2;
  yposy_low  = ypos[1-1][cy] - smallUnit;
  yposy_high = ypos[rows-1][cy] + smallUnit;
  double dx = (yposx_high - yposx_low)/(dimx-1);
  double dy = (yposy_high - yposy_low)/(dimx-1);
  double px, py;

  itk::Array<double> optypos;
  itk::Array<double> deltapos;
  optypos = ypos.get_row(rows-1);
  std::cout << "DebugOptYPos = " << optypos << ";" << std::endl;
  std::cout << "DebugCx = " << cx+1 << ";" << std::endl;
  std::cout << "DebugCy = " << cy+1 << ";" << std::endl;

  for (UnInt i=0; i<dimx; i++)
    {
    px = yposx_low + i * dx;
    std::cout << "DebugYPos1(" << i+1 << ") = " << px << ";" << std::endl;
    }
  for (UnInt j=0; j<dimy; j++)
    {
    py = yposy_low + j * dy;
    std::cout << "DebugYPos2(" << j+1 << ") = " << py << ";" << std::endl;
    }
  for (UnInt i=0; i<dimx; i++)
    {
    px = yposx_low + i * dx;
    for (UnInt j=0; j<dimy; j++)
      {
      py = yposy_low + j * dy;
      ypostmp = optypos;
      ypostmp[cx] = px;
      ypostmp[cy] = py;
      xpostmp = ypostmp * coeffBack + xmean;

      deltapos = xpostmp - this->m_Metric->GetParameters();
      this->m_Metric->UpdateTransformParameters(deltapos);
      this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );

      //std::cout << "DebugPosition(:," << j+1 << "," << i+1 << ") = " << this->m_Metric->GetParameters() << "';" << std::endl;
      std::cout << "DebugValue(" << j+1 << "," << i+1 << ") = " << this->GetValue() << ";" << std::endl;
      }
    }
}

/** L-BFGS Algorithm */
void QuasiNewtonObjectOptimizer
::EstimateNewtonStepWithLBFGS()
{
  const int npar = this->m_Metric->GetNumberOfParameters();
  const int k = m_PastValues.size() - 1;
  if (m_PastValues.size() < npar+1)
    {
    for (unsigned int i=0; i<npar; i++)
      {
      this->m_NewtonStep[i] = this->m_Gradient[i] / this->GetScales()[i];
      }
    return;
    }
  HessianType Hk0;
  Hk0.SetSize(npar, npar);
  ParametersType si(npar);
  ParametersType yi(npar);

  // Initialize Hessian to identity matrix
  Hk0.Fill(0.0f);
  for (int i=0; i<npar; i++)
    {
    Hk0[i][i] = -1.0; //identity matrix //maximizing
    }

  DerivativeType q = m_PastGradients[k];
  std::vector<double> rou(m_PastValues.size());
  std::vector<double> alpha(m_PastValues.size());
  std::vector<double> beta(m_PastValues.size());

  for (int i=k-1; i>=k-npar; i--)
    {
    si = m_PastPositions[i+1] - m_PastPositions[i];
    yi = m_PastGradients[i+1] - m_PastGradients[i];
    rou[i] = inner_product(yi, si);
    alpha[i] = rou[i] * inner_product(si, q);
    q -= alpha[i] * yi;
    }
  ParametersType z(npar);
  z = Hk0 * q;
  for (int i=k-npar; i<=k-1; i++)
    {
    si = m_PastPositions[i+1] - m_PastPositions[i];
    yi = m_PastGradients[i+1] - m_PastGradients[i];
    beta[i] = rou[i] * inner_product(yi, z);
    z += si * (alpha[i] - beta[i]);
    }
  //z = -z; //if minimizing
  this->m_NewtonStep = -z;
}

/** Conjugate gradient method, Fletcher-Reeves Algorithm */
void QuasiNewtonObjectOptimizer
::AdvanceWithConjugateGradient()
{
  const UnInt npar = this->m_Metric->GetNumberOfParameters();
  //x: original parameters, y: scaled parameters
  ParametersType ypos0(npar), ypos1(npar);
  ParametersType xstep(npar), ystep(npar);
  DerivativeType yder0(npar), yder1(npar), ydir(npar);

  for (int loop = 0; loop < 30; loop++)
    {
    if (loop == 0)
      {
      this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
      if (this->GetDebug())
        {
        int iter = 1 + this->GetCurrentIteration();
        std::cout << "m_CurrentPosition(" << iter << ",:) = " << this->m_Metric->GetParameters() << "';" << std::endl;
        std::cout << "m_Value(" << iter << ") = " << this->GetValue() << ";" << std::endl;
        std::cout << "m_Gradient(" << iter << ",:) = " << m_Gradient << "';" << std::endl;
        }
      this->m_CurrentIteration++;

      m_PastPositions.push_back(this->m_Metric->GetParameters());
      m_PastValues.push_back(this->m_Value);
      m_PastGradients.push_back(this->m_Gradient);
      //ypos0 = this->ScalePosition(this->m_Metric->GetParameters());
      this->ScaleDerivative(this->m_Gradient, yder0); // divided by sqrt(scales)
      yder0 = -yder0; //m_Gradient is for maximization, here we want yder for minimization
      ydir = -yder0; //ydir towards optimal
      }
    else
      {
      ydir = -yder0; //ydir towards optimal
      //this->EstimateNewtonStepWithLBFGS();
      //this->ScalePosition(m_NewtonStep, ystep); //multiplied by sqrt(scales)
      //ydir = ystep;
      }

    for (UnInt k=0; k<npar; k++)
      {
      double alpha;

      this->ScaleBackPosition(ydir, xstep); // divided by sqrt(scales)
      alpha = this->EstimateLearningRate( xstep ); //gradient step
      this->AdvanceWithStrongWolfeLineSearch(xstep, alpha*10, alpha);

      //if (loop > 0 && k == 0)
      //  {
      //  alpha = vnl_math_min(alpha, 1.0); //newton step <= 1
      //  }
      xstep *= alpha;

      this->m_Metric->UpdateTransformParameters( xstep );

      this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
      if (this->GetDebug())
        {
        int iter = 1 + this->GetCurrentIteration();
        std::cout << "m_CurrentPosition(" << iter << ",:) = " << this->m_Metric->GetParameters() << "';" << std::endl;
        std::cout << "m_Value(" << iter << ") = " << this->GetValue() << ";" << std::endl;
        std::cout << "m_Gradient(" << iter << ",:) = " << m_Gradient << "';" << std::endl;
        std::cout << "alpha(" << iter << ") = " << alpha << ";" << std::endl;
        std::cout << "xstep(" << iter << ") = " << xstep << ";" << std::endl;
        }
      this->m_CurrentIteration++;

      m_PastPositions.push_back(this->m_Metric->GetParameters());
      m_PastValues.push_back(this->m_Value);
      m_PastGradients.push_back(this->m_Gradient);

      //ypos1 = this->ScalePosition(this->m_Metric->GetParameters());
      this->ScaleDerivative(this->m_Gradient, yder1);
      yder1 = -yder1; //m_Gradient is for maximization, here we want yder for minimization

      double norm1 = yder1.two_norm();
      double norm0 = yder0.two_norm();
      if (norm0 <= NumericTraits<double>::epsilon())
        {
        m_StopCondition = STEP_TOO_SMALL;
        m_StopConditionDescription << "Optimization stops after "
                                   << this->GetCurrentIteration()
                                   << " iterations since"
                                   << " the gradient is too small in line search.";
        this->StopOptimization();
        return;
        }
      double beta = norm1 * norm1 / (norm0 * norm0);
      ydir = -yder1 + beta * ydir;
      if (this->GetDebug())
        {
        int iter = this->GetCurrentIteration();
        std::cout << "norm0(" << iter << ") = " << norm0 << ";" << std::endl;
        std::cout << "beta(" << iter << ") = " << beta << ";" << std::endl;
        }
      yder0 = yder1;
      }
    }
}

/** Do backtracking line search on the Newton direction */
void QuasiNewtonObjectOptimizer
::AdvanceWithBacktrackingLineSearch(ParametersType direction, double maxStepSize,
                                    double &learningRate)
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
      m_StopCondition = STEP_TOO_SMALL;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small in line search.";
      learningRate = t;
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
      m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
      learningRate = t;
      this->StopOptimization();
      return;
      }

    t *= beta;

    }

  //std::cout << "line search step size = " << t << std::endl;
  learningRate = t;
  this->InvokeEvent( IterationEvent() );

}

/** Do line search on a direction with strong Wolfe's conditions.
 * The line search algorithm is from "Introduction to Nonlinear Optimization"
 * by Paul J Atzberger, on
 * Page 7 on http://www.math.ucsb.edu/~atzberg/finance/nonlinearOpt.pdf
 */
void QuasiNewtonObjectOptimizer
::AdvanceWithStrongWolfeLineSearch(ParametersType direction, double maxStepSize,
                                   double &learningRate)
{
  const unsigned int spaceDimension =  this->m_Metric->GetNumberOfParameters();
  double optimalDirection = -1.0; //maximizing

  double tmax = maxStepSize, t0 = 0;
  double t1 = t0, t2 = tmax;

  double c1 = 1e-4, c2 = 0.9;

  double f0, f1, f2;
  double g0, g1, g2; //derivative w.r.t the step size t

  f0 = optimalDirection * this->GetValue();
  g0 = optimalDirection * inner_product(direction, this->m_Gradient);

  f1 = f0;
  g1 = g0;

  int loop = 0;
  ParametersType initPosition = this->m_Metric->GetParameters();
  ParametersType tempPosition(spaceDimension);
  ParametersType deltaPosition(spaceDimension);

  tempPosition = this->m_Metric->GetParameters();
  deltaPosition = initPosition + t2 * direction - tempPosition;
  this->m_Metric->UpdateTransformParameters( deltaPosition );
  this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
  m_CurrentIteration++;

  f2 = optimalDirection * this->m_Value;
  g2 = optimalDirection * inner_product(direction, this->m_Gradient);

  if (-g2 >= -c2 * g0)
    {
    //it may be difficult to find an intermediate t2 such that
    //the stronger wolfe condition holds: vcl_abs(g2) <= -c2 * g0
    if (f2 < f0 + c1 * t2 * g0)
      {
      learningRate = t2;
      return;
      }
    else
      {
      this->AdvanceWithBacktrackingLineSearch(direction, maxStepSize, learningRate);
      return;
      }
    }

  t2 = tmax / 2.0;
  while (true)
    {
    loop++;
    tempPosition = this->m_Metric->GetParameters();
    deltaPosition = initPosition + t2 * direction - tempPosition;

    this->m_Metric->UpdateTransformParameters( deltaPosition );
    this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
    m_CurrentIteration++;

    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = STEP_TOO_SMALL;
      m_StopConditionDescription << "Optimization stops after "
                                 << this->GetCurrentIteration()
                                 << " iterations since"
                                 << " the gradient is too small in line search.";
      this->StopOptimization();

      learningRate = t2;
      return;
      }

    f2 = optimalDirection * this->m_Value;
    g2 = optimalDirection * inner_product(direction, this->m_Gradient);
    if (this->GetDebug())
      {
      std::cout << "LineSearch: phi(" << loop << ",:)=[" << t2 << " " << f2 << " " << g2 << "]" << std::endl;
      }

    if (f2 > f0 + c1 * t2 * g0 || ( loop > 1 && f2 >= f1 ))
      {
      learningRate = this->LineSearchZoom(initPosition, f0, g0, direction, t1, t2);
      return;
      }

    if (vcl_abs(g2) <= -c2 * g0)
      {
      learningRate = t2;
      return;
      }

    if (g2 >= 0)
      {
      learningRate = this->LineSearchZoom(initPosition, f0, g0, direction, t2, t1);
      this->InvokeEvent( IterationEvent() );
      return;
      }
    t1 = t2;
    f1 = f2;
    g1 = g2;

    t2 = t2 + (tmax - t2) / 2;

    /* Update and check iteration count */
    if ( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                                 << ") exceeded in line search. ";
      m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
      this->StopOptimization();
      learningRate = t2;
      return;
      }
    } //while

}

/** Do zoom search on a direction with strong Wolfe's conditions.
 * The zoom algorithm is from "Introduction to Nonlinear Optimization"
 * by Paul J Atzberger, on
 * Page 7 on http://www.math.ucsb.edu/~atzberg/finance/nonlinearOpt.pdf
 */
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
      std::cout << "LineSearchZoom: m_CurrentPosition=" << initPosition + t2 * direction << std::endl;
      std::cout << "LineSearchZoom: m_Value=" << this->m_Value << std::endl;
      std::cout << "LineSearchZoom: m_Gradient=" << this->m_Gradient << std::endl;
      }

    double gradientNorm = m_Gradient.two_norm();
    if (gradientNorm < m_MinimumGradientNorm)
      {
      m_StopCondition = STEP_TOO_SMALL;
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
      m_StopCondition = STEP_TOO_SMALL;
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
      m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
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

/** Estimate the quasi-newton step */
void QuasiNewtonObjectOptimizer
::EstimateNewtonStep()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  // Estimate Hessian
  EstimateHessian();

  // Compute the Newton step
  ParametersType sdx(numPara);

  if (this->GetCurrentIteration() == 0 ||
      //this->GetCurrentIteration() == 100 ||
      m_Hessian[0][0] == NumericTraits<double>::max() ||
      m_HessianInverse[0][0] == NumericTraits<double>::max())
    {
    this->ResetNewtonStep();
    }
  else
    {
    sdx = m_HessianInverse * m_Gradient;
    //this->m_NewtonStep = -1.0 * sdx; //minimize
    this->m_NewtonStep = sdx; //maximize
    }

  if ( this->GetDebug() )
    {
    unsigned int iter = this->GetCurrentIteration()+1;
    for (unsigned int row = 0; row < m_Hessian.rows(); row++)
      {
      std::cout << "m_Hessian(" << row+1 << ",:" << "," << iter << ") = ["
        << m_Hessian.get_row(row) << "];" << std::endl;
      }
    }

}

/** Estimate Hessian matrix with BFGS method described
 *  at http://en.wikipedia.org/wiki/BFGS_method
 */
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
    //HessianType approximateHessian(numPara, numPara);
    //approximateHessian.set_identity();
    //approximateHessian *= m_Lambda;
    //approximateHessian += newHessian;
    //m_HessianInverse = vnl_matrix_inverse<double>(approximateHessian);

    m_HessianInverse = vnl_matrix_inverse<double>(newHessian);
    }
}

void QuasiNewtonObjectOptimizer
::ResetNewtonStep()
{
  unsigned int numPara = this->m_Metric->GetNumberOfParameters();

  //m_NewtonStep = m_Gradient; //use gradient step
  for (unsigned int i=0; i<numPara; i++)
    {
    this->m_NewtonStep[i] = this->m_Gradient[i] / this->GetScales()[i];
    }

  // Initialize Hessian to identity matrix
  m_Hessian.Fill(0.0f);
  m_HessianInverse.Fill(0.0f);

  for (unsigned int i=0; i<numPara; i++)
    {
    m_Hessian[i][i] = 1.0; //identity matrix
    m_HessianInverse[i][i] = 1.0; //identity matrix
    }
}

/** Compute the learning late from voxel shift*/
double QuasiNewtonObjectOptimizer
::EstimateLearningRate(ParametersType step, double maxshift)
{
  if (m_OptimizerParameterEstimator.IsNull())
    {
    return 1;
    }

  ParametersType parameters = this->m_CurrentPosition;

  ScalesType     scales = this->GetScales();

  double shift, learningRate;

  shift = m_OptimizerParameterEstimator->ComputeMaximumVoxelShift(step);

  if (this->GetCurrentIteration() == 0)
    {
    std::cout << " Initial learningRate = " << m_MaximumVoxelShift / shift << std::endl;
    }

  if (shift > 0)
    {
    //learningRate = m_MaximumVoxelShift / shift;
    learningRate = maxshift / shift;
    }
  else //if the step yields almost zero voxel shift, it is not good to go
    {
    learningRate = 0;
    }
  return learningRate;
}
/** Translate the parameters into the scaled space */
void QuasiNewtonObjectOptimizer
::ScalePosition(const ParametersType p1, ParametersType &p2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<p1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      p2[i] = p1[i] * scale;
    }
}

/** Translate the parameters into the original space */
void QuasiNewtonObjectOptimizer
::ScaleBackPosition(const ParametersType p1, ParametersType &p2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<p1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      p2[i] = p1[i] / scale;
    }
}

/** Translate the derivative into the scaled space */
void QuasiNewtonObjectOptimizer
::ScaleDerivative(const DerivativeType g1, DerivativeType &g2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<g1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      g2[i] = g1[i] / scale;
    }
}

/** Translate the derivative into the original space */
void QuasiNewtonObjectOptimizer
::ScaleBackDerivative(const DerivativeType g1, DerivativeType &g2)
{
  double scale = 1.0;
  for (unsigned int i=0; i<g1.size(); i++)
    {
      scale = vcl_sqrt(this->GetScales()[i]);
      g2[i] = g1[i] * scale;
    }
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
