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
//#ifndef __itkQuasiNewtonLocalSupportObjectOptimizer_hxx
//#define __itkQuasiNewtonLocalSupportObjectOptimizer_hxx

#include <iostream>
#include <iomanip>

#include "itkQuasiNewtonLocalSupportObjectOptimizer.h"
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
QuasiNewtonLocalSupportObjectOptimizer
::QuasiNewtonLocalSupportObjectOptimizer()
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
QuasiNewtonLocalSupportObjectOptimizer
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * Start the optimization
 */
void
QuasiNewtonLocalSupportObjectOptimizer
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
        itkExceptionMacro("QuasiNewtonLocalSupportObjectOptimizer: metric parameters are not defined.");
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
QuasiNewtonLocalSupportObjectOptimizer
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
      m_StopCondition = COSTFUNCTION_ERROR;
      m_StopConditionDescription << "CostFunction error";
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
QuasiNewtonLocalSupportObjectOptimizer
::AdvanceOneStep(void)
{
  itkDebugMacro("AdvanceOneStep");

  if ( this->m_Metric->HasLocalSupport() )
    {
    AdvanceOneLocalStep();
    return;
    }
}

/**
 * Advance one Step: if the Quasi-Newton direction is consistent
 * with the gradient direction, follow the Quasi-Newton direction,
 * otherwise, follow the gradient direction.
 */
void
QuasiNewtonLocalSupportObjectOptimizer
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
void QuasiNewtonLocalSupportObjectOptimizer
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
void QuasiNewtonLocalSupportObjectOptimizer
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
double QuasiNewtonLocalSupportObjectOptimizer
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
void QuasiNewtonLocalSupportObjectOptimizer
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
