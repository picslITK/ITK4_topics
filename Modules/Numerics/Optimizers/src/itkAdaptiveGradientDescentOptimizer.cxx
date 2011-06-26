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
#ifndef _itkAdaptiveGradientDescentOptimizer_txx
#define _itkAdaptiveGradientDescentOptimizer_txx

#include "itkSigmoidImageFilter.h"
#include "itkAdaptiveGradientDescentOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkMacro.h"

namespace itk
{
/**
 * Constructor
 */
AdaptiveGradientDescentOptimizer
::AdaptiveGradientDescentOptimizer()
{
  itkDebugMacro("Constructor");

  m_Param_a = 1;
  m_Param_A = 20;
  m_Param_alpha = 0.6;

  m_CurrentTime = 0.0;
  m_InitialTime = 0.0;

  m_SigmoidMax = 1.0;
  m_SigmoidMin = -0.8;
  m_SigmoidScale = 1e-8;

}

void
AdaptiveGradientDescentOptimizer
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "OptimizerHelper: "
     << m_OptimizerHelper << std::endl;
}

/**
 * Advance one Step following the gradient direction
 */
void
AdaptiveGradientDescentOptimizer
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

  ScalesType scales(this->GetScales().size());
  double shiftLearningRate;
  double adaptiveLearningRate;

  // initialize the parameter scales
  if (this->GetCurrentIteration() == 0)
    {
    m_PreviousGradient = m_Gradient;
    }

  adaptiveLearningRate = this->EstimateLearningRate( );

  if (this->GetCurrentIteration() == 0)
    {
    // adjust the m_Param_a such that the first step will produce a shift
    // of at most MaximumVoxelShift
    ParametersType step = m_Gradient;
    step = direction * step;
    shiftLearningRate = this->Superclass::EstimateLearningRate(step);
    m_Param_a = m_Param_a * shiftLearningRate / adaptiveLearningRate;
    adaptiveLearningRate = this->EstimateLearningRate( );
    }

  this->SetLearningRate( adaptiveLearningRate );

  this->GradientDescentOptimizer::AdvanceOneStep();

  this->UpdateCurrentTime();

}


/**
* This function computes the parameter a at iteration, as
* described by Spall.
*/
double AdaptiveGradientDescentOptimizer
::EstimateLearningRate() const
{
  double base = m_Param_A + m_CurrentTime + 1.0;

  return m_Param_a / vcl_pow( base, m_Param_alpha );

}

/**
 * Update the m_CurrentTime by m_Gradient and m_PreviousGradient
 * using a sigmoid function.
 */
void AdaptiveGradientDescentOptimizer
  ::UpdateCurrentTime( void )
{
  typedef itk::Functor::Sigmoid<double, double> SigmoidType;

  if ( this->GetCurrentIteration() > 0 )
    {
    SigmoidType sigmoid;
    sigmoid.SetOutputMaximum( m_SigmoidMax );
    sigmoid.SetOutputMinimum( m_SigmoidMin );
    sigmoid.SetAlpha( m_SigmoidScale );

    double beta = m_SigmoidScale *
       vcl_log( - m_SigmoidMax / m_SigmoidMin );
    sigmoid.SetBeta( beta );

    double prod = inner_product(
      this->m_PreviousGradient, this->GetGradient() );

    this->m_CurrentTime = this->m_CurrentTime + sigmoid(-prod);
    this->m_CurrentTime = vnl_math_max( 0.0, this->m_CurrentTime );
    }

  this->m_PreviousGradient = this->GetGradient();

} // end UpdateCurrentTime


} // end namespace itk

#endif
