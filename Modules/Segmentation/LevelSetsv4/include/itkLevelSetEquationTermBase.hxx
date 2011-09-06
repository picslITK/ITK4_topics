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

#ifndef __itkLevelSetEquationTermBase_hxx
#define __itkLevelSetEquationTermBase_hxx

#include "itkLevelSetEquationTermBase.h"

#include "itkNumericTraits.h"

namespace itk
{
// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::LevelSetEquationTermBase()
{
  this->m_CurrentLevelSetId =  NumericTraits< LevelSetIdentifierType >::Zero;
  this->m_Coefficient = NumericTraits< LevelSetOutputRealType >::One;
  this->m_CFLContribution = NumericTraits< LevelSetOutputRealType >::Zero;
  this->m_TermName = "";
}

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::~LevelSetEquationTermBase()
{
}

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
const typename LevelSetEquationTermBase< TInputImage, TLevelSetContainer >::RequiredDataType &
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::GetRequiredData() const
{
  return this->m_RequiredData;
}

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
void
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::SetLevelSetContainer( LevelSetContainerType* iContainer )
{
  this->m_LevelSetContainer = iContainer;
  this->m_Heaviside = iContainer->GetHeaviside();
  this->Modified();
}

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
typename
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::LevelSetOutputRealType
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::Evaluate( const LevelSetInputIndexType& iP )
{
  return this->m_Coefficient * this->Value( iP );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
typename
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::LevelSetOutputRealType
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::Evaluate( const LevelSetInputIndexType& iP,
            const LevelSetDataType& iData )
{
  return this->m_Coefficient * this->Value( iP, iData );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TInputImage, class TLevelSetContainer >
void
LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
::SetUp()
{
  this->m_CFLContribution = NumericTraits< LevelSetOutputRealType >::Zero;

  if( this->m_CurrentLevelSetPointer.IsNull() )
    {
    this->m_CurrentLevelSetPointer =
    this->m_LevelSetContainer->GetLevelSet( this->m_CurrentLevelSetId );

    if( this->m_CurrentLevelSetPointer.IsNull() )
      {
      itkWarningMacro(
      << "m_CurrentLevelSetId does not exist in the level set container" );
      }
    }

  if( !this->m_Heaviside.IsNotNull() )
    {
    itkWarningMacro( << "m_Heaviside is NULL" );
    }
}
// ----------------------------------------------------------------------------

}

#endif // __itkLevelSetEquationTermBase_hxx
