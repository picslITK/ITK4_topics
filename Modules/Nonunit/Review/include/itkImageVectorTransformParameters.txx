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

#ifndef __itkImageVectorTransformParameters_txx
#define __itkImageVectorTransformParameters_txx

#include "itkImageVectorTransformParameters.h"

namespace itk
{
/** Default contstructor */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::ImageVectorTransformParameters() : TransformParameters< TValueType >()
{
  m_ParameterImage = NULL;
}

/** Copy constructor */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::ImageVectorTransformParameters(const ImageVectorTransformParameters& rhs)
  : TransformParameters< TValueType >(rhs)
{
  m_ParameterImage = NULL;
}

/** Constructor with size */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::ImageVectorTransformParameters(unsigned int dimension)
  : TransformParameters< TValueType >(dimension)
{
  m_ParameterImage = NULL;
}

/** Move the data pointer */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
void
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::MoveDataPointer( TValueType * pointer )
{
  // The buffer for Image<Vector> points to Vector type, not TValueType, so
  // have to cast.
  typedef typename ParameterImageType::PixelContainer::Element vectorElement;
  vectorElement* vectorPointer = reinterpret_cast<vectorElement *>(pointer);
  // We're expecting the new memory buffer t be of same size.
  unsigned int sizeInVectors = m_ParameterImage->GetPixelContainer()->Size();
  // After this call, PixelContainer will *not* manage its memory.
  this->m_ParameterImage->GetPixelContainer()->SetImportPointer( vectorPointer,
                                                              sizeInVectors );
  Superclass::MoveDataPointer( pointer );
}

/** Set parameter image */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
void
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::SetParameterImage( ParameterImageType * image)
{
  m_ParameterImage = image;
  //The PixelContainer for Image<Vector> points to type Vector, so we have
  // to determine the number of raw elements of type TValueType in the buffer
  // and cast a pointer to it for assignment to the Array data pointer.
  unsigned int sz = image->GetPixelContainer()->Size() * NVectorDimension;
  TValueType* valuePointer = reinterpret_cast<TValueType *>
                            ( image->GetPixelContainer()->GetBufferPointer() );
  //Set the Array's pointer to the image data buffer
  this->SetData( valuePointer, sz );
}

template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
const typename ImageVectorTransformParameters< TValueType,
                                               NVectorDimension,
                                               VImageDimension >
::Self &
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::operator=(const Self & rhs)
{
  if ( this == &rhs ) { return *this; }

  // Call the superclass implementation
  this->Superclass::operator=(rhs);

  return *this;
}

template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
const typename ImageVectorTransformParameters< TValueType,
                                               NVectorDimension,
                                               VImageDimension >
::Self &
ImageVectorTransformParameters< TValueType, NVectorDimension, VImageDimension >
::operator=(const Superclass & rhs)
{
  if ( this == &rhs ) { return *this; }

  // Call the superclass implementation
  this->Superclass::operator=(rhs);

  return *this;
}


}//namespace itk
#endif
