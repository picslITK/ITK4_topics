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
  this->m_ValueAndDerivativeThreader->SetHolder( this );

  /* These will be instantiated if needed in Initialize */
  m_MovingGradientCalculator = NULL;
  m_FixedGradientCalculator = NULL;
}

/*
 * Initialize
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::Initialize()
{
  //Allocate and initialize memory for holding derivative results.
  this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
  this->m_MeasurePerThread.resize( this->m_NumberOfThreads );
  /* This size always comes from the moving image */
  unsigned long globalDerivativeSize =
    this->m_Metric->GetMovingImageTransform()->GetNumberOfParameters();
  itkDebugMacro("ImageToImageObjectMetric::Initialize: deriv size  "
                  << globalDerivativeSize << std::endl);
  this->m_GlobalDerivative.SetSize( globalDerivativeSize );
  /* For transforms with local support, e.g. deformation field,
   * use a single derivative container that's updated by region
   * in multiple threads. */
  if ( this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx HAS local support\n");
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetData(
                                    this->m_GlobalDerivative.data_block(),
                                    this->m_GlobalDerivative.Size(),
                                    false );
      }
    }
  else
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx does NOT have local support\n");
    /* Global transforms get a separate derivatives container for each thread
     * that holds the result for a particular region. */
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetSize( globalDerivativeSize );
      /* Be sure to init to 0 here, because the threader may not use
       * all the threads if the region is better split into fewer
       * subregions. */
      this->m_DerivativesPerThread[i].Fill( 0 );
      }
    }
}

/*
 * Threader callback.
 * Iterates over image region and calls derived class for calculations.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeMultiThreadedCallback(
                const ThreaderInputObjectType& virtualImageSubRegion,
                int threadId,
                Self * self)
{
  /* Create an iterator over the virtual sub region */
  ImageRegionConstIteratorWithIndex<VirtualImageType>
      ItV( this->m_VirtualDomainImage, virtualImageSubRegion );

  VirtualPointType            virtualPoint;
  FixedOutputPointType        mappedFixedPoint;
  FixedImagePixelType         fixedImageValue;
  FixedImageDerivativesType   fixedImageDerivatives;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        movingImageValue;
  MovingImageDerivativesType  movingImageDerivatives;
  bool                        pointIsValid;
  MeasureType                 metricValue;

  //FIXME: handle different coord systems, will need diff tranform?
  //  what if using coord system 'both'?
  DerivativeType            localDerivative( self->GetLocalDerivativeSize() );

  /* Iterate over the sub region */
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {

    self->m_VirtualImage->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);

    /* Transform the point into fixed and moving spaces, and evaluate.
     * These methods will check that the point lies within the mask if
     * one has been set, and then verify they lie in the fixed or moving
     * space as appropriate.
     * If both tests pass, the point is evaluated and 'true' is returned. */
    self->TransformAndEvaluateFixedPoint( virtualPoint,
                                          mappedFixedPoint,
                                          pointIsValid,
                                          fixedImageValue,
                                          true /*compute gradient*/,
                                          fixedImageDerivatives,
                                          threadID );
    if( pointIsValid )
      {
      self->TransformAndEvaluateMovingPoint( virtualPoint,
                                            mappedMovingPoint,
                                            pointIsValid,
                                            movingImageValue,
                                            true /*compute gradient*/,
                                            movingImageDerivatives,
                                            threadID );
      }

    /* Call the user method in derived classes to do the specific
     * calculations for value and derivative. */
    if( pointIsValid )
      {
      pointIsValid = self->GetValueAndDerivativeProcessPoint(
                                               virtualPoint,
                                               mappedFixedPoint,
                                               fixedImageValue,
                                               fixedImageDerivatives,
                                               mappedMovingPoint,
                                               movingImageValue,
                                               movingImageDerivatives,
                                               metricValue,
                                               localDerivative,
                                               threadID );
      }

    if( pointIsValid )
      {
      increment count of pixels used in computation, prob in per-thread var
      }

    //next index
    ++ItV;
  }

  //save results to per-thread containers
}

/*
 * Transform a point from VirtualImage domain to FixedImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint( VirtualPointType & point,
                      FixedImagePointType & mappedFixedPoint,
                      bool & pointIsValid,
                      FixedImagePixelType & fixedImageValue,
                      bool computeImageGradient,
                      FixedImageDerivativesType & fixedGradient,
                      unsigned int threadID) const
{
TODO: make sure image deriv is done right in terms of doing in either fixed or moving space and then mapping to virutal, like the switch I made in Demons. Need to multiply by Jacobian after transforming into virtual space, Brian mentioned this. Am I doing that already?

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
        if( computeImageGradient )
          {
          /* The assumption here is that calling EvaluateValueAndDerivative
           * is faster than separately calling Evaluate and EvaluateDerivative,
           * because of setup for the bspline evaluation. */
          m_FixedBSplineInterpolator->EvaluateValueAndDerivative(
                                                        mappedFixedPoint,
                                                        fixedImageValue,
                                                        fixedGradient,
                                                        threadID);
          }
        else
          {
          fixedImageValue =
            m_FixedBSplineInterpolator->Evaluate(mappedFixedPoint, threadID);
          }
        }
      }
    else
      {
      // Check if mapped point inside image buffer
      pointIsValid = m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
      if ( pointIsValid )
        {
        fixedImageValue = m_FixedInterpolator->Evaluate(mappedFixedPoint);
        if( computeImageGradient )
          {
          this->ComputeFixedImageDerivatives(mappedFixedPoint,
                                             fixedGradient,
                                             threadID);
          }
        }
      }
    }
}

/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint( VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      bool computeImageGradient,
                      MovingImageDerivativesType & movingGradient,
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
        if( computeImageGradient )
          {
          /* The assumption here is that calling EvaluateValueAndDerivative
           * is faster than separately calling Evaluate and EvaluateDerivative,
           * because of setup for the bspline evaluation. */
          m_MovingBSplineInterpolator->EvaluateValueAndDerivative(
                                                        mappedMovingPoint,
                                                        movingImageValue,
                                                        movingGradient,
                                                        threadID);
          }
        else
          {
          movingImageValue =
            m_MovingBSplineInterpolator->Evaluate(mappedMovingPoint, threadID);
          }
        }
      }
    else
      {
      // Check if mapped point inside image buffer
      pointIsValid = m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
      if ( pointIsValid )
        {
        movingImageValue = m_MovingInterpolator->Evaluate(mappedMovingPoint);
        if( computeImageGradient )
          {
          this->ComputeMovingImageDerivatives(mappedMovingPoint,
                                              movingGradient,
                                              threadID);
          }
        }
      }
    }
}

/*
 * Compute image derivatives for a Fixed point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageDerivatives(
                              const FixedImagePointType & mappedPoint,
                              FixedImageDerivativesType & gradient,
                              unsigned int threadID) const
{
  if ( m_FixedInterpolatorIsBSpline )
    {
    // Computed Fixed image gradient using derivative BSpline kernel.
    gradient = m_FixedBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                              threadID);
    }
  else
    {
    if ( m_UseGaussianGradient )
      {
      ContinuousIndex< double, FixedImageDimension > tempIndex;
      m_FixedImage->TransformPhysicalPointToContinuousIndex(mappedPoint,
                                                             tempIndex);
      FixedImageIndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);
      gradient = m_FixedGaussianGradientImage->GetPixel(mappedIndex);
      }
    else
      {
      // if not using the gradient image
      gradient = m_FixedGradientCalculator->Evaluate(mappedPoint);
      }
    }
}

/*
 * Compute image derivatives for a moving point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageDerivatives(
                              const MovingImagePointType & mappedPoint,
                              MovingImageDerivativesType & gradient,
                              unsigned int threadID) const
{
  if ( m_MovingInterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_MovingBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                               threadID);
    }
  else
    {
    if ( m_UseGaussianGradient )
      {
      ContinuousIndex< double, MovingImageDimension > tempIndex;
      m_MovingImage->TransformPhysicalPointToContinuousIndex(mappedPoint,
                                                             tempIndex);
      MovingImageIndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);
      gradient = m_MovingGaussianGradientImage->GetPixel(mappedIndex);
      }
    else
      {
      // if not using the gradient image
      gradient = m_MovingGradientCalculator->Evaluate(mappedPoint);
      }
    }
}

/*
 * Get the size of localDerivative
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
SizeValueType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetLocalDerivativeSize()
{
  switch( this->m_CoordinateSystem )
    {
    case Moving:
      return this->m_MovingImageTransform->GetNumberOfLocalParameters() );
      break;
    case Fixed:
      return this->m_FixedImageTransform->GetNumberOfLocalParameters() );
      break;
    case Both:
      /* TODO: this isn't completely settled, what to do here. Derived
       * classes may need to override to provide different behavior. */
      return VirtualImageDimension;
      break;
    default:
      itkExceptionMacro("Invalid CoordinateSystem enum: "
                         << this->m_CoordinateSystem );
    }
  return 0;
}

}//namespace itk
