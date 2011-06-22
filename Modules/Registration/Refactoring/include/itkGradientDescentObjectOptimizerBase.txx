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
::ModifyGradientThreaded( const IndexRangeType& rangeForThread,
                          int threadId,
                          Self *holder )
{
  holder->ModifyGradientOverSubRange( rangeForThread );
}

}//namespace itk

#endif
