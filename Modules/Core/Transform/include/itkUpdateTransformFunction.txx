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
#ifndef __itkUpdateTransformFunction_txx
#define __itkUpdateTransformFunction_txx

#include "itkUpdateTransformFunction.h"

namespace itk
{

template <class TTransform>
void
UpdateTransformFunction<TTransform>
::Update( DerivativeType & update,
          ScalarType factor,
          TransformType * transform )
{
  unsigned int numberOfParameters = transform->GetNumberOfParameters();
  if( update.Size() != numberOfParameters )
    {
    itkExceptionMacro("Parameter update size, " << update.Size() << ", must "
                      " be same as transform parameter size, "
                      << numberOfParameters << std::endl);
    }
  if( factor == 1.0 )
    {
    for (unsigned int k=0; k < numberOfParameters; k++)
      transform->m_Parameters[k] += update[k];
    }
  else
    {
    for (unsigned int k=0; k < numberOfParameters; k++)
      transform->m_Parameters[k] += update[k] * factor;
    }
}

}

#endif
