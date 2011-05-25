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
#ifndef __itkGradientDescentThreadedMetricOptimizer_txx
#define __itkGradientDescentThreadedMetricOptimizer_txx

#include "itkGradientDescentThreadedMetricOptimizer.h"

namespace itk
{
/** Advance one step in the optimization */
template<class TMetricFunction, class TThreader>
void
GradientDescentThreadedMetricOptimizer<TMetricFunction,TThreader>
::AdvanceOneStep()
{
  //ScalesType &    scales = this->GetScales();

}


}//namespace itk

#endif
