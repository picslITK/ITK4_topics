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
#ifndef __itkImageToSpatialObjectMetric_txx
#define __itkImageToSpatialObjectMetric_txx

#include "itkImageToSpatialObjectMetric.h"

namespace itk
{
/** Constructor */
template< class TFixedImage, class TMovingSpatialObject, typename TValueType >
ImageToSpatialObjectMetric< TFixedImage, TMovingSpatialObject, TValueType >
::ImageToSpatialObjectMetric()
{
  m_FixedImage          = 0; // has to be provided by the user.
  m_MovingSpatialObject = 0; // has to be provided by the user.
  m_Transform           = 0; // has to be provided by the user.
  m_Interpolator        = 0; // has to be provided by the user.
}

/**
 * Initialize
 */

template< class TFixedImage, class TMovingSpatialObject, typename TValueType >
void
ImageToSpatialObjectMetric< TFixedImage, TMovingSpatialObject, TValueType >
::Initialize(void)
throw ( ExceptionObject )
{
  if ( !m_Transform )
    {
    itkExceptionMacro(<< "Transform is not present");
    }

  if ( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator is not present");
    }

  if ( !m_MovingSpatialObject )
    {
    itkExceptionMacro(<< "MovingSpatialObject is not present");
    }

  if ( !m_FixedImage )
    {
    itkExceptionMacro(<< "FixedImage is not present");
    }

  // If the image is provided by a source, update the source.
  if ( m_FixedImage->GetSource() )
    {
    m_FixedImage->GetSource()->Update();
    }

  m_Interpolator->SetInputImage(m_FixedImage);

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
}

/** PrintSelf */
template< class TFixedImage, class TMovingSpatialObject, typename TValueType >
void
ImageToSpatialObjectMetric< TFixedImage, TMovingSpatialObject, TValueType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving Spatial Object: " << m_MovingSpatialObject.GetPointer()  << std::endl;
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer()   << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Last Transform parameters = " << m_LastTransformParameters << std::endl;
}
} // end namespace itk

#endif
