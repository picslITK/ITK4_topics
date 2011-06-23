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
#ifndef __itkGradientDescentObjectOptimizer_txx
#define __itkGradientDescentObjectOptimizer_txx

#include "itkGradientDescentObjectOptimizer.h"

namespace itk
{

/**
 * Default constructor
 */
void
GradientDescentObjectOptimizer
::GradientDescentObjectOptimizer()
{
  m_LearningRate = 1.0;
}

/**
 * Start and run the optimization
 */
void
GradientDescentObjectOptimizer
::StartOptimization()
{
  itkDebugMacro("StartOptimization");

  m_CurrentIteration   = 0;

  /* Validate some settings */
  if( this->m_Metric.IsNull() )
    {
    itkExceptionMacro("m_Metric must be set.");
    return;
    }

  if( this->m_Scales.Size() != this->m_Metric->GetNumberOfParameters() )
    {
    m_StopCondition = OtherError;
    m_StopConditionDescription << "Scales size error";
    StopOptimization();
    itkExceptionMacro("Size of scales (" << this->m_Scales.Size()
                      << ") must match size of parameters (" <<
                      this->m_Metric->GetNumberOfParameters() << ").");
    }

  this->ResumeOptimization();
}

/**
 * Resume optimization.
 */
void
GradientDescentObjectOptimizer
::ResumeOptimization()
{
  /* Do threading initialization here so we can also do some cleanup
   * when we stop, and still be able to resume. */
  this->InitializeForThreading();

  m_StopConditionDescription.str("");
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
  InvokeEvent( StartEvent() );

  m_Stop = false;
  while( ! m_Stop )
    {
    /* Compute value/derivative, using threader. */
    try
      {
      this->m_Metric->GetValueAndDerivative( this->m_Value, this->m_Gradient );
      }
    catch ( ExceptionObject & err )
      {
      m_StopCondition = MetricError;
      m_StopConditionDescription << "Metric error";
      StopOptimization();

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
      StopOptimization();
      break;
      }
    } //while (!m_Stop)
}

/**
 * Advance one Step following the gradient direction
 */
void
GradientDescentObjectOptimizer
::AdvanceOneStep()
{
  itkDebugMacro("AdvanceOneStep");

  /* Begin threaded gradient modification. Work is done in
   * ModifyGradientOverSubRange */
  this->ModifyGradient();

  this->m_Metric->UpdateTransformParameters( m_Gradient );

  this->InvokeEvent( IterationEvent() );
}

/**
 * Modify the gradient over a given index range.
 */
void
GradientDescentObjectOptimizer
::ModifyGradientOverSubRange( IndexRangeType& subrange )
{
  double direction;
  if ( this->m_Maximize )
    {
    direction = 1.0;
    }
  else
    {
    direction = -1.0;
    }

  ScalesType& scales = this->GetScales();

  double direction_learning = direction * m_LearningRate;
  for ( unsigned int j = subrange[0]; j <= subrange[1]; j++ )
    {
    m_Gradient[j] = m_Gradient[j] / scales[j] * direction_learning;
    }
}

}//namespace itk

#endif
