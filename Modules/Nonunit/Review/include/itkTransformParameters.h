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
#ifndef __itkTransformParameters_h
#define __itkTransformParameters_h

#include "itkArray.h"

namespace itk
{
/** \class TransformParameters
 *  \brief Class to hold and manage different parameter types used by Transforms.
 *
 */

template< typename TValueType >
class TransformParameters : public Array< TValueType >
{
public:

  /** The element type stored at each location in the Array. */
  typedef TValueType                          ValueType;
  typedef TransformParameters                 Self;
  typedef Array< TValueType >                 Superclass;
  typedef Superclass                          ArrayType;
  typedef typename Superclass::VnlVectorType  VnlVectorType;

  /** Default constructor. It is created with an empty array
   *  it has to be allocated later by assignment              */
  TransformParameters();

  /** Copy constructor.  Uses VNL copy construtor with correct
   *  setting for memory management.
   *  The vnl vector copy constructor creates new memory
   *  no matter the setting of let array manage memory of rhs.
   */
  TransformParameters(const TransformParameters& rhs);

  /** Constructor with size. Size can only be changed by assignment */
  explicit TransformParameters(unsigned int dimension);

  /** Set a new data pointer for the Array, pointing it to a different
   * memory block.
   * The size of the new memroy block must be the same,
   * in elements of TValueType.
   * Memory must be managed by caller afterwards. */
  virtual void MoveDataPointer( TValueType * pointer );

  /** Copy opertor */
  const Self & operator=(const Self & rhs);

  const Self & operator=(const ArrayType & rhs);

  const Self & operator=(const VnlVectorType & rhs);

  ~TransformParameters(){}
};

}//namespace itk

#if ITK_TEMPLATE_TXX
#include "itkTransformParameters.txx"
#endif

#endif
