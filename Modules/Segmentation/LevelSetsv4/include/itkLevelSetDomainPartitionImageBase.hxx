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
#ifndef __itkLevelSetDomainPartitionImageBase_hxx
#define __itkLevelSetDomainPartitionImageBase_hxx

#include "itkLevelSetDomainPartitionImageBase.h"

namespace itk
{
template< class TImage >
LevelSetDomainPartitionImageBase< TImage >
::LevelSetDomainPartitionImageBase()
{
}

template< class TImage >
LevelSetDomainPartitionImageBase< TImage >
::~LevelSetDomainPartitionImageBase()
{
}

template< class TImage >
void LevelSetDomainPartitionImageBase< TImage >
::PopulateListDomain()
{
  const ListRegionType & region = this->m_ListDomain->GetLargestPossibleRegion();
  ListIteratorType lIt(this->m_ListDomain, region);

  for ( lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
    {
    ListIndexType listIndex = lIt.GetIndex();
    IdentifierListType identifierList;
    IdentifierType i = NumericTraits< IdentifierType >::Zero;
    while( i < this->m_NumberOfLevelSetFunctions )
      {
      if ( this->m_LevelSetDataPointerVector[i]->VerifyInsideRegion(listIndex) )
        {
        identifierList.push_back(i);
        }
      i++;
      }
    lIt.Set(identifierList);
    }
}

template< class TImage >
void LevelSetDomainPartitionImageBase< TImage >
::AllocateListDomain()
{
  if( m_Image.IsNull() )
    {
    itkGenericExceptionMacro( "m_Image is null" );
    }
  this->m_ListDomain = ListImageType::New();
  this->m_ListDomain->CopyInformation( this->m_Image );
  this->m_ListDomain->SetRegions( this->m_Image->GetLargestPossibleRegion() );
  this->m_ListDomain->Allocate();
}

} //end namespace itk

#endif
