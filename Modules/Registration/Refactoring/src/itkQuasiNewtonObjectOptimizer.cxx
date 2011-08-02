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

  // m_MinimumVoxelShift will be set to a value in the first iteration
  // according the first step
  m_MinimumVoxelShift = 0;

  m_LocalHessian = NULL;
  m_LocalHessianInverse = NULL;
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

  if ( ! this->m_Metric->HasLocalSupport() )
    {
    // initialize scales
    ScalesType scales(this->m_Metric->GetNumberOfParameters());
    //this->m_Metric->EstimateScales(true, scales);
    m_OptimizerParameterEstimator->EstimateScales(this->m_Metric->GetParameters(), scales);
    this->SetScales(scales);
    std::cout << " Estimated scales = " << scales << std::endl;
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

  //std::cout << "m_NewtonStep = " << m_NewtonStep << std::endl;
  //std::cout << "m_Gradient = "   << m_Gradient << std::endl;

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
    //std::cout << "using gradient, learningRate = " << learningRate << std::endl;

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
    //std::cout << "using newton, learningRate = " << learningRate << std::endl;

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

  //this->SetDebug(true);
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
    m_NewtonStep.Fill(0); //use gradient step

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
        m_NewtonStep[i*imageDimension + d] = 0; //use gradient
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
  ParametersType parameters = this->GetCurrentPosition();

  ScalesType     scales = this->GetScales();

  double shift, learningRate;

  //shift = this->m_Metric->ComputeMaximumVoxelShift(true, step);
  shift = m_OptimizerParameterEstimator->ComputeMaximumVoxelShift(parameters, step);

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

//#endif
