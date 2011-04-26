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

#ifndef __itkTransformParameters_txx
#define __itkTransformParameters_txx

#include "itkTransformParameters.h"

namespace itk
{
/** Default contstructor */
template< typename TValueType >
TransformParameters< TValueType >
::TransformParameters() : Array< TValueType >()
{
}

/** Copy constructor */
template< typename TValueType >
TransformParameters< TValueType >
::TransformParameters(const TransformParameters& rhs)
  : Array< TValueType >(rhs)
{
}

/** Constructor with size */
template< typename TValueType >
TransformParameters< TValueType >
::TransformParameters(unsigned int dimension)
  : Array< TValueType >(dimension)
{
}

template< typename TValueType >
void
TransformParameters< TValueType >
::MoveDataPointer( TValueType * pointer)
{
  this->SetData( pointer, this->GetSize(), false /*LetArrayManageMemory*/);
}

template< typename TValueType >
const typename TransformParameters< TValueType >
::Self &
TransformParameters< TValueType >
::operator=(const Self & rhs)
{
  if ( this == &rhs ) { return *this; }

  // Call the superclass implementation
  this->ArrayType::operator=(rhs);

  return *this;
}

template< typename TValueType >
const typename TransformParameters< TValueType >
::Self &
TransformParameters< TValueType >
::operator=(const ArrayType & rhs)
{
  if ( this == &rhs ) { return *this; }

  // Call the superclass implementation
  this->ArrayType::operator=(rhs);

  return *this;
}

template< typename TValueType >
const typename TransformParameters< TValueType >
::Self &
TransformParameters< TValueType >
::operator=(const VnlVectorType & rhs)
{
  if ( this == &rhs ) { return *this; }

  // Call the superclass implementation
  this->ArrayType::operator=(rhs);

  return *this;
}

}//namespace itk
#endif
