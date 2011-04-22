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
#ifndef __itkMembershipFunctionBase_h
#define __itkMembershipFunctionBase_h

#include "itkFunctionBase.h"
#include "itkMeasurementVectorTraits.h"
#include "itkNumericTraitsCovariantVectorPixel.h"

namespace itk
{
namespace Statistics
{
/** \class MembershipFunctionBase
 * \brief MembershipFunctionBase class declares common interfaces
 * for membership functions.
 *
 * As a function derived from FunctionBase, users use Evaluate method
 * get result. However, the return value type of the method is fixed
 * as double. Any function derived from this class returns quantitative
 * measure for how well the vector x belong to the class ( or group)
 * represented by the function.
 * \ingroup ITK-Statistics
 */

template< class TVector >
class ITK_EXPORT MembershipFunctionBase:
  public FunctionBase< TVector, double >
{
public:
  /** Standard class typedefs */
  typedef MembershipFunctionBase          Self;
  typedef FunctionBase< TVector, double > Superclass;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;

  /** Strandard macros */
  itkTypeMacro(MembershipFunctionBase, FunctionBase);

  /** MeasurementVector typedef support */
  typedef TVector MeasurementVectorType;

  /** Typedef for the length of each measurement vector */
  typedef unsigned int MeasurementVectorSizeType;

  /** Method to get membership score (discriminant score) of an entity. */
  virtual double Evaluate(const MeasurementVectorType & x) const = 0;

  /** Set method for the length of the measurement vector */
  virtual void SetMeasurementVectorSize(MeasurementVectorSizeType s)
  {
    // Test whether the vector type is resizable or not
    MeasurementVectorType m;

    if ( MeasurementVectorTraits::IsResizable(m) )
      {
      // then this is a resizable vector type
      //
      // if the new size is the same as the previou size, just return
      if ( s == this->m_MeasurementVectorSize )
        {
        return;
        }
      else
        {
        this->m_MeasurementVectorSize = s;
        this->Modified();
        }
      }
    else
      {
      // If this is a non-resizable vector type
      MeasurementVectorType     m3;
      MeasurementVectorSizeType defaultLength =
        NumericTraits<MeasurementVectorType>::GetLength(m3);
      // and the new length is different from the default one, then throw an
      // exception
      if ( defaultLength != s )
        {
        itkExceptionMacro(
          "Attempting to change the measurement \
           vector size of a non-resizable vector type" );
        }
      }
  }

  /** Get method for the length of the measurement vector */
  itkGetConstMacro(MeasurementVectorSize, MeasurementVectorSizeType);
protected:
  MembershipFunctionBase()
  {
    m_MeasurementVectorSize = NumericTraits<MeasurementVectorType>::GetLength(
      MeasurementVectorType() );
  }

  virtual ~MembershipFunctionBase(void) {}

  void PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "Length of measurement vectors: "
       << m_MeasurementVectorSize << std::endl;
  }

  MeasurementVectorSizeType m_MeasurementVectorSize;
};  // end of class
} // end of namespace Statistics
} // end namespace itk

#endif
