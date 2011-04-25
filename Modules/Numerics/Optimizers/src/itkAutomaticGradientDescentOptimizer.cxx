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
#ifndef _itkAutomaticGradientDescentOptimizer_txx
#define _itkAutomaticGradientDescentOptimizer_txx

#include "itkAutomaticGradientDescentOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkMacro.h"

namespace itk
{
/**
 * Constructor
 */
AutomaticGradientDescentOptimizer
::AutomaticGradientDescentOptimizer()
{
  itkDebugMacro("Constructor");

  m_MaximumVoxelShift = 1.0;
}

void
AutomaticGradientDescentOptimizer
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "m_MaximumVoxelShift: "
     << m_MaximumVoxelShift << std::endl;
  os << indent << "OptimizerHelper: "
     << m_OptimizerHelper << std::endl;
}

/**
 * Start the optimization
 */
void
AutomaticGradientDescentOptimizer
::StartOptimization(void)
{
  itkDebugMacro("StartOptimization");

  // initialize scales
  ScalesType scales(this->GetInitialPosition().Size());
  m_OptimizerHelper->EstimateScales(this->GetInitialPosition(), scales);
  this->SetScales(scales);
  std::cout << " Estimated scales = " << scales << std::endl;

  // m_MinimumVoxelShift will be set to a value in the first iteration
  // according the first step
  m_MinimumVoxelShift = 0;

  this->Superclass::StartOptimization();

}

/**
 * Advance one Step following the gradient direction
 */
void
AutomaticGradientDescentOptimizer
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

  // initialize learningRate
  if (this->GetCurrentIteration() == 0)
    {
    double learningRate;
    ParametersType step = m_Gradient;
    step = direction * step;
    learningRate = this->EstimateLearningRate(step);
    this->SetLearningRate( learningRate );
    }

  this->Superclass::AdvanceOneStep();

}

/** Compute learning late from voxel shift*/
double AutomaticGradientDescentOptimizer
::EstimateLearningRate(ParametersType step)
{
  ParametersType parameters = this->GetCurrentPosition();
  int            numPara    = parameters.size();

  ScalesType     scales = this->GetScales();

  ParametersType deltaParameters(numPara);
  for (int i=0; i<numPara; i++)
    {
      deltaParameters[i] = step[i] / scales[i];
    }

  double shift, learningRate;

  shift = m_OptimizerHelper->ComputeMaximumVoxelShift(parameters, deltaParameters);

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
void AutomaticGradientDescentOptimizer
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
