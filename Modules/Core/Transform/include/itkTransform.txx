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
#ifndef __itkTransform_txx
#define __itkTransform_txx

#include "itkTransform.h"

namespace itk
{
/**
 * Constructor
 */
template< class TScalarType,
          unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
Transform< TScalarType, NInputDimensions, NOutputDimensions >
::Transform():
  m_Parameters(1),
  m_FixedParameters(1),
  m_Jacobian(NOutputDimensions, 1)
{
  m_IdentityJacobian.SetSize(NOutputDimensions,NOutputDimensions);
  m_IdentityJacobian.Fill(0);
  for( unsigned int i=0; i < NOutputDimensions; i++ )
    {
    m_IdentityJacobian[i][i] = 1.0;
    }

  /* Set the default parameters update function */
  m_UpdateTransformFunction = UpdateTransformFunctionType::New();

  itkWarningMacro(
    << "Using default transform constructor.  Should specify NOutputDims and NParameters as args to constructor.");
}

/**
 * Constructor
 */
template< class TScalarType,
          unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
Transform< TScalarType, NInputDimensions, NOutputDimensions >
::Transform(unsigned int dimension, unsigned int numberOfParameters):
  m_Parameters(numberOfParameters),
  m_FixedParameters(numberOfParameters),
  m_Jacobian(dimension, numberOfParameters)
{
  m_IdentityJacobian.SetSize(NOutputDimensions,NOutputDimensions);
  m_IdentityJacobian.Fill(0);
  for( unsigned int i=0; i < NOutputDimensions; i++ )
    {
    m_IdentityJacobian[i][i] = 1.0;
    }
  /* Set the default parameters update function */
  m_UpdateTransformFunction = UpdateTransformFunctionType::New();
}

/**
 * GenerateName
 */
template< class TScalarType,
          unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
std::string Transform< TScalarType, NInputDimensions, NOutputDimensions >
::GetTransformTypeAsString() const
{
  std::ostringstream n;

  n << GetNameOfClass();
  n << "_";
  n << this->GetTransformTypeAsString(static_cast<TScalarType *>(0));
  n << "_" << this->GetInputSpaceDimension() << "_" << this->GetOutputSpaceDimension();
  return n.str();
}

/**
 * UpdateTransformParameters
 */
template< class TScalarType,
          unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
Transform< TScalarType, NInputDimensions, NOutputDimensions >
::UpdateTransformParameters( DerivativeType & update,
                              TScalarType factor )
{
  /* The default function adds the update after scaling by factor. */
  this->m_UpdateTransformFunction->Update( update, factor, this );

  /* Call SetParameters with the updated parameters.
   * SetParameters in most transforms is used to assign the input params
   * to member variables, possibly with some processing. The member variables
   * are then used in TransformPoint.
   * In the case of dense-field transforms that are updated in blocks from
   * a threaded implementation, SetParameters doesn't do this, and is
   * optimized to not copy the input parameters when == m_Parameters.
   */
  this->SetParameters( this->m_Parameters );

  /* Call Modified, following behavior of other transform when their
   * parameters change, e.g. MatrixOffsetTransformBase */
  this->Modified();
}

} // end namespace itk

#endif
