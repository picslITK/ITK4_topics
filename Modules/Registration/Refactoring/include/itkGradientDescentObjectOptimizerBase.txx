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
#ifndef __itkGradientDescentObjectOptimizerBase_txx
#define __itkGradientDescentObjectOptimizerBase_txx

#include "itkGradientDescentObjectOptimizerBase.h"

namespace itk
{

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::GradientDescentObjectOptimizerBase()
{
  this->m_ModifyGradientThreader->SetThreadedGenerateData(
    Self::ModifyGradientThreaded );

  m_Maximize = false;
  m_NumberOfIterations = 100;
  m_CurrentIteration = 0;
  m_StopCondition = MaximumNumberOfIterations;
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
}

//-------------------------------------------------------------------
template<class TMetricFunction>
const std::string
GradientDescentObjectOptimizerBase<TMetricFunction>
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
const std::string
GradientDescentObjectOptimizerBase<TMetricFunction>
::StopOptimization(void)
{
  itkDebugMacro("StopOptimization");

  m_Stop = true;
  this->CleanupFromThreading();
  InvokeEvent( EndEvent() );
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::InitializeForThreading()
{
  if( this->m_Metric.IsNull() )
    {
    itkExceptionMacro("m_Metric must be set.");
    return;
    }

  // TODO: call m_Metric->Initialize(), or check that
  // it's been init'ed somehow? Would all metrics have Initialize(), or
  // at least all new ones to be used with this optimzier hierarchy?

/*

for regular step grad descent:

  const unsigned int dimensions =
    this->m_Metric->GetNumberOfParameters();
  // Verify parameter settings
  if ( m_RelaxationFactor < 0.0 )
    {
    itkExceptionMacro(<< "Relaxation factor must be positive. "
                         "Current value is " << m_RelaxationFactor);
    return;
    }

  if ( m_RelaxationFactor >= 1.0 )
    {
    itkExceptionMacro(<< "Relaxation factor must less than 1.0. "
                         "Current value is " << m_RelaxationFactor);
    return;
    }

*/
  std::cout << " end Initialize " << std::endl;
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::UpdateMetricValueAndDerivative()
{
//  this->BeforeMetricThreadedGenerateData();
  /* Calculate new metric value and derivative with threading. */
  this->m_MetricThreader->GenerateData();
  /* Collect metric results from the threads. */
//  this->AfterMetricThreadedGenerateData();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::BeforeMetricThreadedGenerateData()
{
  /* Allow the metric to do any per-iteration initialization it
   * requires. This method is not expected to be thread-safe, and
   * may implement its own threading */
//  this->m_Metric->BeforeThreadedGetValueAndDerivative();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::AfterMetricThreadedGenerateData()
{

//  std::cout << " end after threaded generate data. global derivative: "
  //          << this->m_GlobalDerivative << std::endl;
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::CleanupFromThreading()
{
removed metric-threading stuff
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentObjectOptimizerBase<TMetricFunction>
::ModifyGradientThreaded( const IndexRangeType& rangeForThread,
                          int threadId,
                          Self *holder )
{
  holder->ModifyGradientOverSubRange( rangeForThread );
}

}//namespace itk

#endif
