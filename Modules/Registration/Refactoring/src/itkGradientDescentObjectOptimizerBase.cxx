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
//#ifndef __itkGradientDescentObjectOptimizerBase_txx
//#define __itkGradientDescentObjectOptimizerBase_txx

#include "itkGradientDescentObjectOptimizerBase.h"

namespace itk
{

//-------------------------------------------------------------------
GradientDescentObjectOptimizerBase
::GradientDescentObjectOptimizerBase()
{
  this->m_ModifyGradientThreader = ModifyGradientThreaderType::New();
  this->m_ModifyGradientThreader->SetThreadedGenerateData(
    Self::ModifyGradientThreaded );
  this->m_ModifyGradientThreader->SetHolder( this );

  m_NumberOfIterations = 100;
  m_CurrentIteration = 0;
  m_StopCondition = MaximumNumberOfIterations;
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
}

//-------------------------------------------------------------------
const std::string
GradientDescentObjectOptimizerBase
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

//-------------------------------------------------------------------
void
GradientDescentObjectOptimizerBase
::StopOptimization(void)
{
  itkDebugMacro("StopOptimization");
  m_Stop = true;
  InvokeEvent( EndEvent() );
}

//-------------------------------------------------------------------
void
GradientDescentObjectOptimizerBase
::ModifyGradient()
{
  IndexRangeType fullrange;
  fullrange[0] = 0;
  fullrange[1] = this->m_Gradient.GetSize()-1; //range is inclusive
  /* Perform the modification either with or without threading */
  if( this->m_Metric->HasLocalSupport() )
    {
    this->m_ModifyGradientThreader->SetOverallIndexRange( fullrange );
    /* This ends up calling ModifyGradientThreaded from each thread */
    this->m_ModifyGradientThreader->GenerateData();
    }
  else
    {
    /* Global transforms are small, so update without threading. */
    this->ModifyGradientOverSubRange( fullrange );
    }
}
//-------------------------------------------------------------------
void
GradientDescentObjectOptimizerBase
::ModifyGradientThreaded( const IndexRangeType& rangeForThread,
                          ThreadIdType threadId,
                          Self *holder )
{
  holder->ModifyGradientOverSubRange( rangeForThread );
}

}//namespace itk

//#endif
