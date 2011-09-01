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
#ifndef __itkImageAdaptor_hxx
#define __itkImageAdaptor_hxx
#include "itkImageAdaptor.h"
#include "itkProcessObject.h"

namespace itk
{
/**
 * Constructor
 */
template< class TImage, class TAccessor >
ImageAdaptor< TImage, TAccessor >
::ImageAdaptor()
{
  // Allocate an internal image.  A process object might try to allocate an
  // temporary image that is the same type as its input or output.  If that
  // image type is an adaptor, we need to make sure that an internal image is
  // available because the process object will not know to call SetImage on
  // the adaptor.
  m_Image = TImage::New();
}

/**
 * Destructor
 */
template< class TImage, class TAccessor >
ImageAdaptor< TImage, TAccessor >
::~ImageAdaptor()
{}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::Initialize()
{
  // call the superclass' method first; then delegate
  Superclass::Initialize();

  // delegation to internal image
  m_Image->Initialize();
}

template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetPixelContainer(PixelContainer *container)
{
  if ( this->GetPixelContainer() != container )
    {
    m_Image->SetPixelContainer(container);
    this->Modified();
    }
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::Graft(const DataObject *data)
{
  // call the superclass' implementation
  Superclass::Graft(data);

  if ( data )
    {
    // Attempt to cast data to an ImageAdaptor
    const Self *imgData;

    try
      {
      imgData = dynamic_cast< const Self * >( data );
      }
    catch ( ... )
      {
      return;
      }

    if ( imgData )
      {
      // Now copy anything remaining that is needed
      this->SetPixelContainer(
        const_cast< Self * >( imgData )->GetPixelContainer() );
      }
    else
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::ImageAdaptor::Graft() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( const Self * ).name() );
      }
    }
}

/**
 *
 */
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::OffsetValueType *
ImageAdaptor< TImage, TAccessor >
::GetOffsetTable(void) const
{
  return m_Image->GetOffsetTable();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
typename ImageAdaptor< TImage, TAccessor >::IndexType
ImageAdaptor< TImage, TAccessor >
::ComputeIndex(OffsetValueType offset) const
{
  return m_Image->ComputeIndex(offset);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::Update()
{
  Superclass::Update();

  m_Image->Update();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::UpdateOutputInformation()
{
  // call the superclass' method first, then delegate
  Superclass::UpdateOutputInformation();

  // delegation to internal image
  m_Image->UpdateOutputInformation();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::UpdateOutputData()
{
  // call the superclass' method first, then delegate
  Superclass::UpdateOutputData();

  // delegation to internal image
  m_Image->UpdateOutputData();
  SetBufferedRegion( m_Image->GetBufferedRegion() );
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::PropagateRequestedRegion()
throw ( InvalidRequestedRegionError )
{
  // call the superclass' method first, then delegate
  Superclass::PropagateRequestedRegion();

  // delegation to internal image
  m_Image->PropagateRequestedRegion();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetRequestedRegionToLargestPossibleRegion()
{
  // call the superclass' method first, then delegate
  Superclass::SetRequestedRegionToLargestPossibleRegion();

  // delegation to internal image
  m_Image->SetRequestedRegionToLargestPossibleRegion();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::CopyInformation(const DataObject *data)
{
  // call the superclass' method first, then delegate
  Superclass::CopyInformation(data);

  // delegation to internal image
  m_Image->CopyInformation(data);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::SpacingType &
ImageAdaptor< TImage, TAccessor >
::GetSpacing(void) const
{
  return m_Image->GetSpacing();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetSpacing(const SpacingType & spacing)
{
  // delegation to internal image
  m_Image->SetSpacing(spacing);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetSpacing(const double *spacing /*[Self::ImageDimension]*/)
{
  // delegation to internal image
  m_Image->SetSpacing(spacing);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetSpacing(const float *spacing /*[Self::ImageDimension]*/)
{
  // delegation to internal image
  m_Image->SetSpacing(spacing);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetOrigin(const PointType origin)
{
  // delegation to internal image
  m_Image->SetOrigin(origin);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetOrigin(const double *origin /*[Self::ImageDimension]*/)
{
  // delegation to internal image
  m_Image->SetOrigin(origin);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetOrigin(const float *origin /*[Self::ImageDimension]*/)
{
  // delegation to internal image
  m_Image->SetOrigin(origin);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::PointType &
ImageAdaptor< TImage, TAccessor >
::GetOrigin(void) const
{
  return m_Image->GetOrigin();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetDirection(const DirectionType direction)
{
  // delegation to internal image
  m_Image->SetDirection(direction);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::DirectionType &
ImageAdaptor< TImage, TAccessor >
::GetDirection(void) const
{
  return m_Image->GetDirection();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetImage(TImage *image)
{
  m_Image = image;
  Superclass::SetLargestPossibleRegion( m_Image->GetLargestPossibleRegion() );
  Superclass::SetBufferedRegion( m_Image->GetBufferedRegion() );
  Superclass::SetRequestedRegion( m_Image->GetRequestedRegion() );
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::InternalPixelType  *
ImageAdaptor< TImage, TAccessor >
::GetBufferPointer() const
{
  return m_Image->GetBufferPointer();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
typename ImageAdaptor< TImage, TAccessor >::InternalPixelType  *
ImageAdaptor< TImage, TAccessor >
::GetBufferPointer()
{
  return m_Image->GetBufferPointer();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::Modified() const
{
  Superclass::Modified();

  m_Image->Modified();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
unsigned long
ImageAdaptor< TImage, TAccessor >
::GetMTime() const
{
  unsigned long mtime1, mtime2;

  mtime1 = Superclass::GetMTime();
  mtime2 = m_Image->GetMTime();

  return ( mtime1 >= mtime2 ? mtime1 : mtime2 );
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetBufferedRegion(const RegionType & region)
{
  // call the superclass' method first, then delegate
  Superclass::SetBufferedRegion(region);

  // delegation to internal image
  m_Image->SetBufferedRegion(region);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::RegionType &
ImageAdaptor< TImage, TAccessor >
::GetBufferedRegion(void) const
{
  // delegation to internal image
  return m_Image->GetBufferedRegion();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetLargestPossibleRegion(const RegionType & region)
{
  // call the superclass' method first, then delegate
  Superclass::SetLargestPossibleRegion(region);

  // delegation to internal image
  m_Image->SetLargestPossibleRegion(region);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::RegionType &
ImageAdaptor< TImage, TAccessor >
::GetLargestPossibleRegion(void) const
{
  // delegation to internal image
  return m_Image->GetLargestPossibleRegion();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetRequestedRegion(const RegionType & region)
{
  // call the superclass' method first, then delegate
  Superclass::SetRequestedRegion(region);

  // delegation to internal image
  m_Image->SetRequestedRegion(region);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
void
ImageAdaptor< TImage, TAccessor >
::SetRequestedRegion(const DataObject *data)
{
  // call the superclass' method first, then delegate
  Superclass::SetRequestedRegion(data);

  // delegation to internal image
  m_Image->SetRequestedRegion(data);
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
bool
ImageAdaptor< TImage, TAccessor >
::VerifyRequestedRegion()
{
  // call the superclass' method first, then delegate
  Superclass::VerifyRequestedRegion();

  // delegation to internal image
  return m_Image->VerifyRequestedRegion();
}

//----------------------------------------------------------------------------
template< class TImage, class TAccessor >
const typename ImageAdaptor< TImage, TAccessor >::RegionType &
ImageAdaptor< TImage, TAccessor >
::GetRequestedRegion(void) const
{
  // delegation to internal image
  return m_Image->GetRequestedRegion();
}
} // end namespace itk

#endif
