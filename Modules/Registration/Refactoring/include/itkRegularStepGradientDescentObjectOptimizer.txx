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
#ifndef __itkRegularStepGradientDescentObjectOptimizer_txx
#define __itkRegularStepGradientDescentObjectOptimizer_txx

#include "itkRegularStepGradientDescentObjectOptimizer.h"

namespace itk
{
/** Advance one step in the optimization */
template<class TMetricFunction>
void
RegularStepGradientDescentObjectOptimizer<TMetricFunction>
::AdvanceOneStep()
{
  //ScalesType &    scales = this->GetScales();

}


/*

for regular step grad descent initialize:

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

}//namespace itk

#endif
