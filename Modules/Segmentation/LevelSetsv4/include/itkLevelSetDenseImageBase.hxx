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

#ifndef __itkLevelSetDenseImageBase_hxx
#define __itkLevelSetDenseImageBase_hxx

#include "itkLevelSetDenseImageBase.h"

namespace itk
{
// ----------------------------------------------------------------------------
template< class TImage >
LevelSetDenseImageBase< TImage >
::LevelSetDenseImageBase()
{}

// ----------------------------------------------------------------------------
template< class TImage >
LevelSetDenseImageBase< TImage >
::~LevelSetDenseImageBase()
{
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetDenseImageBase< TImage >
::SetImage( ImageType* iImage )
{
  this->m_Image = iImage;
  typename ImageType::SpacingType spacing = m_Image->GetSpacing();

  for( unsigned int dim = 0; dim < Dimension; dim++ )
    {
    this->m_NeighborhoodScales[dim] =
      NumericTraits< OutputRealType >::One / static_cast< OutputRealType >( spacing[dim ] );
    }
  this->Modified();
}

// ----------------------------------------------------------------------------
template< class TImage >
typename LevelSetDenseImageBase< TImage >::OutputType
LevelSetDenseImageBase< TImage >::Evaluate( const InputType& iP ) const
{
  return this->m_Image->GetPixel( iP );
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetDenseImageBase< TImage >::Evaluate( const InputType& iP, LevelSetDataType& ioData ) const
{
  Superclass::Evaluate( iP, ioData );
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetDenseImageBase< TImage >
::Initialize()
{
  Superclass::Initialize();

  this->m_Image = NULL;
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetDenseImageBase< TImage >
::CopyInformation(const DataObject *data)
{
  Superclass::CopyInformation( data );

  const Self *LevelSet = NULL;

  try
    {
    LevelSet = dynamic_cast< const Self * >( data );
    }
  catch ( ... )
    {
    // LevelSet could not be cast back down
    itkExceptionMacro( << "itk::LevelSetDenseImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }

  if ( !LevelSet )
    {
    // pointer could not be cast back down
    itkExceptionMacro( << "itk::LevelSetDenseImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetDenseImageBase< TImage >
::Graft( const DataObject* data )
{
  Superclass::Graft( data );
  const Self *LevelSet = NULL;

  try
    {
    LevelSet = dynamic_cast< const Self* >( data );
    }
  catch( ... )
    {
    // image could not be cast back down
    itkExceptionMacro( << "itk::LevelSetDenseImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
    }

  if ( !LevelSet )
    {
    // pointer could not be cast back down
    itkExceptionMacro( << "itk::LevelSetDenseImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }

  this->m_Image = LevelSet->m_Image;
  this->m_NeighborhoodScales = LevelSet->m_NeighborhoodScales;
}

// ----------------------------------------------------------------------------
template< class TImage >
bool
LevelSetDenseImageBase< TImage >
::IsInside(const InputType &iP) const
{
  const RegionType largestRegion = this->m_Image->GetLargestPossibleRegion();

  return largestRegion.IsInside( iP );
}

}
#endif // __itkLevelSetDenseImageBase_hxx
