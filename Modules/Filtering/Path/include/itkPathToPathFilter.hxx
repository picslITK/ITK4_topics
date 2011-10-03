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
#ifndef __itkPathToPathFilter_hxx
#define __itkPathToPathFilter_hxx

#include "itkPathToPathFilter.h"

namespace itk
{
/**
 *
 */
template< class TInputPath, class TOutputPath >
PathToPathFilter< TInputPath, TOutputPath >
::PathToPathFilter()
{
  // Let the superclass do everything
}

/**
 *
 */
template< class TInputPath, class TOutputPath >
void
PathToPathFilter< TInputPath, TOutputPath >
::SetInput(const InputPathType *path)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< InputPathType * >( path ) );
}

/**
 * Connect one of the operands for pixel-wise addition
 */
template< class TInputPath, class TOutputPath >
void
PathToPathFilter< TInputPath, TOutputPath >
::SetInput(unsigned int index, const TInputPath *path)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( index,
                                    const_cast< TInputPath * >( path ) );
}

/**
 *
 */
template< class TInputPath, class TOutputPath >
const typename PathToPathFilter< TInputPath, TOutputPath >::InputPathType *
PathToPathFilter< TInputPath, TOutputPath >
::GetInput(void)
{
  return static_cast< const TInputPath * >( this->GetPrimaryInput() );
}

/**
 *
 */
template< class TInputPath, class TOutputPath >
const typename PathToPathFilter< TInputPath, TOutputPath >::InputPathType *
PathToPathFilter< TInputPath, TOutputPath >
::GetInput(unsigned int idx)
{
  return static_cast< const TInputPath * >
         ( this->ProcessObject::GetInput(idx) );
}

/**
 *
 */
template< class TInputPath, class TOutputPath >
void
PathToPathFilter< TInputPath, TOutputPath >
::GenerateInputRequestedRegion()
{
  // ProcessObject::GenerateInputRequestedRegion() will (for each input) call
  // Path::SetRequestedRegionToLargestPossibleRegion(), which is empty.
  Superclass::GenerateInputRequestedRegion();
}

/**
 *
 */
template< class TInputPath, class TOutputPath >
void
PathToPathFilter< TInputPath, TOutputPath >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
