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
#ifndef __itkImageToImageObjectMetric_txx
#define __itkImageToImageObjectMetric_txx

#include "itkImageToImageObjectMetric.h"
//#include "itkImageRandomConstIteratorWithIndex.h"

namespace itk
{

/*
 * Constructor
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ImageToImageObjectMetric
{
  this->m_ValueAndDerivativeThreader->SetThreadedGenerateData(
    Self::GetValueAndDerivativeMultiThreadedCallback );
  this->m_ValueAndDerivativeThreader->SetHolder(
}

/*
 * Threader callback.
 * Iterates over image region and calls derived class for calculations.
 */
void GetValueAndDerivativeMultiThreadedCallback(
                const ThreaderInputObjectType& virtualImageRegion,
                int threadId,
                Self * metric)
{
  //loop over virtual region

  //transform points to fixed and moving

  //if both points valid, pass to derived class for calcs

  //save results to per-thread containers
}

/** Transform a point from VirtualImage domain to FixedImage domain.
 * This function also checks if mapped point is within the mask if
 * on is set, and that is within the moving image buffer.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformFixedPoint(VirtualPointType & point,
                      FixedImagePointType & mappedFixedPoint,
                      bool & pointIsValid,
                      FixedImagePixelType & fixedImageValue,
                      unsigned int threadID) const
{
  pointIsValid = true;
  fixedImageValue = 0;
  mappedFixedPoint = m_FixedImageTransform->TransformPoint( point );

  // If user provided a mask over the Moving image
  if ( m_FixedImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_FixedImageMask->IsInside( mappedFixedPoint );
    }

  if( pointIsValid )
    {
    if ( m_FixedInterpolatorIsBSpline )
      {
      // Check if mapped point inside image buffer
      pointIsValid =
        m_FixedBSplineInterpolator->IsInsideBuffer( mappedFixedPoint );
      if ( pointIsValid )
        {
        fixedImageValue =
          m_FixedBSplineInterpolator->Evaluate(mappedFixedPoint, threadID);
        }
      }
    else
      {
      // Check if mapped point inside image buffer
      pointIsValid = m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
      if ( pointIsValid )
        {
        fixedImageValue = m_FixedInterpolator->Evaluate(mappedFixedPoint);
        }
      }
    }
}

/** Transform a point from VirtualImage domain to MovingImage domain.
 * This function also checks if mapped point is within the mask if
 * on is set, and that is within the moving image buffer.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformMovingPoint(VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      unsigned int threadID) const
{
  pointIsValid = true;
  movingImageValue = 0;
  mappedMovingPoint = m_MovingImageTransform->TransformPoint( point );

  // If user provided a mask over the Moving image
  if ( m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_MovingImageMask->IsInside( mappedMovingPoint );
    }

  if( pointIsValid )
    {
    if ( m_MovingInterpolatorIsBSpline )
      {
      // Check if mapped point inside image buffer
      pointIsValid =
        m_MovingBSplineInterpolator->IsInsideBuffer( mappedMovingPoint );
      if ( pointIsValid )
        {
        movingImageValue =
          m_MovingBSplineInterpolator->Evaluate(mappedMovingPoint, threadID);
        }
      }
    else
      {
      // Check if mapped point inside image buffer
      pointIsValid = m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
      if ( pointIsValid )
        {
        movingImageValue = m_MovingInterpolator->Evaluate(mappedMovingPoint);
        }
      }
    }
}

}//namespace itk
