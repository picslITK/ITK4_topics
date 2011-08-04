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

#ifndef __itkImageVectorTransformParametersHelper_hxx
#define __itkImageVectorTransformParametersHelper_hxx

#include "itkImageVectorTransformParametersHelper.h"

namespace itk
{
/** Default contstructor */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
ImageVectorTransformParametersHelper< TValueType, NVectorDimension, VImageDimension >
::ImageVectorTransformParametersHelper()
{
  m_ParameterImage = NULL;
}

/** Move the data pointer */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
void
ImageVectorTransformParametersHelper< TValueType, NVectorDimension, VImageDimension >
::MoveDataPointer( CommonContainerType* container, TValueType * pointer )
{
  if( m_ParameterImage.IsNull() )
    {
    itkGenericExceptionMacro("ImageVectorTransformParametersHelper::"
      "MoveDataPointer: m_ParameterImage must be defined.");
    }
  // The buffer for Image<Vector> points to Vector type, not TValueType, so
  // have to cast.
  typedef typename ParameterImageType::PixelContainer::Element vectorElement;
  vectorElement* vectorPointer = reinterpret_cast<vectorElement *>(pointer);
  // We're expecting the new memory buffer t be of same size.
  unsigned int sizeInVectors = m_ParameterImage->GetPixelContainer()->Size();
  // After this call, PixelContainer will *not* manage its memory.
  this->m_ParameterImage->GetPixelContainer()->SetImportPointer( vectorPointer,
                                                              sizeInVectors );
  Superclass::MoveDataPointer( container, pointer );
}

/** Set parameter image */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
void
ImageVectorTransformParametersHelper< TValueType, NVectorDimension, VImageDimension >
::SetParametersObject(CommonContainerType * container, LightObject * object)
{
  if( object == NULL )
    {
    m_ParameterImage = NULL;
    return;
    }
  else
    {
    ParameterImageType* image =
      dynamic_cast<ParameterImageType *>( object );
    if( image == NULL )
      {
      itkGenericExceptionMacro(
        "ImageVectorTransformParametersHelper::SetParametersObject: object is "
        "not of proper image type. Expected VectorImage, received "
        << object->GetNameOfClass() )
      }
    m_ParameterImage = image;
    //The PixelContainer for Image<Vector> points to type Vector, so we have
    // to determine the number of raw elements of type TValueType in the buffer
    // and cast a pointer to it for assignment to the Array data pointer.
    unsigned int sz = image->GetPixelContainer()->Size() * NVectorDimension;
    TValueType* valuePointer = reinterpret_cast<TValueType *>
                              ( image->GetPixelContainer()->GetBufferPointer() );
    //Set the Array's pointer to the image data buffer. By default it will
    // not manage the memory.
    container->SetData( valuePointer, sz );
    }
}

}//namespace itk
#endif
